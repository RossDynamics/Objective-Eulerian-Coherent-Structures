function Compute_OECSv1
clear all

%Global Variables
global x y u v Eigenvalue1 Eigenvalue2 Eigenvector1 Eigenvector2 index Eigencheck xmin xmax ymin ymax

%Input Parameters
xmin=-2;
xmax=2;
ymin=-1;
ymax=1;
res=10;
dist=1/res;
xlength=((xmax-xmin)*res+1);
ylength=((ymax-ymin)*res+1);

%Calculate Velocity Field & Coordinate Mesh
[x,y,u,v]=velocity(xmin,xmax,ymin,ymax,res);

%Calculate Gradients of Velocity Field
[dxu,dyu]=gradient(u,dist);
[dxv,dyv]=gradient(v,dist);

%Compute Eulerian Rate of Strain Tensor for each point in the Coordinate Mesh
for i=[1:1:ylength]
    for j=[1:1:xlength]
        S(i,j,1,1)=1/2*(dxu(i,j)+dxu(i,j));
        S(i,j,1,2)=1/2*(dxv(i,j)+dyu(i,j));
        S(i,j,2,1)=1/2*(dyu(i,j)+dxv(i,j));
        S(i,j,2,2)=1/2*(dyv(i,j)+dyv(i,j));
    end
end

%Calculate Eigenvalue and Eigenvector fields for Eulerian Strain Tensor
for i=[1:1:ylength]
    for j=[1:1:xlength]
        Matrix(1,1)=S(i,j,1,1);
        Matrix(1,2)=S(i,j,1,2);
        Matrix(2,1)=S(i,j,2,1);
        Matrix(2,2)=S(i,j,2,2);
        A=eig(Matrix);
        [V D]=eig(Matrix);
        if A(1)<A(2)
            Eigenvalue1(i,j)=A(1);
            Eigenvalue2(i,j)=A(2);
            Eigenvector1(i,j,:)=V(:,1);
            Eigenvector2(i,j,:)=V(:,2);
        else
            Eigenvalue1(i,j)=A(2);
            Eigenvalue2(i,j)=A(1);
            Eigenvector1(i,j,:)=V(:,2);
            Eigenvector2(i,j,:)=V(:,1);
        end
        clear Matrix A V D
    end
end
%figure
%hold on
%surf(x,y,Eigenvalue(:,:,2))


%{
%Find Singularities in the Rate of Strain Tensor
%[Xq,Yq] = meshgrid(0:0.01:2,0:0.01:1);
%S1q = interp2(x,y,S(:,:,1,1)-S(:,:,2,2),Xq,Yq,'linear');
%S2q = interp2(x,y,S(:,:,1,2),Xq,Yq,'linear');
%for j=[1:1:201]
%    for i=[1:1:101]

% Find and Classify Singularities
for i=[1:1:ylength]
    for j=[1:1:xlength]
        %if point (j,i) is Singular
        if and(and(S(i,j,1,1)-S(i,j,2,2)<1e-5,S(i,j,1,2)<1e-5),abs(Eigenvalue(i,j,1)*Eigenvalue(i,j,2)-1)>1)
            SingularityMatrix(i,j)=1;
            %ZZ(i,j)=ClassifySingularity([Eigenvector(i,j,1,1),Eigenvector(i,j,1,2)]);
            %XX(i,j)=ClassifySingularity([Eigenvector(i,j,2,1),Eigenvector(i,j,2,2)]);
        else
            SingularityMatrix(i,j)=0;
            %ZZ(i,j)=-1;
            %XX(i,j)=-1;
        end
    end
end

%SingularityMatrix
%figure
%surf(Xq,Yq,SingularityMatrix)
figure
surf(x,y,SingularityMatrix)
%surf(x,y,ZZ)
%figure
%surf(x,y,XX)
end
%}

%Compute Hyperbolic OECS
%Set up grid vectors
x1=x(1,:);
y1=y(:,1);
[s1xmin s1ymin]=Find2DPeak(abs(Eigenvalue1(:,:)),x1,y1,'maxima');
%figure
%hold on
%
%scatter(s1xmin, s1ymin,'r','filled')

[s2xmax s2ymax]=Find2DPeak(abs(Eigenvalue2(:,:)),x1,y1,'maxima');
%figure
%hold on
%quiver(x,y,u,v)


%

Ax0=s2xmax(1);
Ay0=s2ymax(1);
A0=[Ax0,Ay0,0];
Rx0=s1xmin(1);
Ry0=s1ymin(1);
R0=[Rx0,Ry0,0];
timeA=[0 200];
timeR=[0 200];
%EigcheckA
%EigcheckA
%opts=odeset('RelTol',1e-9,'Events', @eventfun);

%opts=odeset('MaxStep',0.1,'Events', @eventfun);
opts=odeset('Events', @eventfun);
clear Eigencheck;
index=0;
[T A1]=ode45(@Attracting1,timeA,A0,opts);
clear Eigencheck;
index=0;
[T A2]=ode45(@Attracting2,timeA,A0,opts);
clear Eigencheck;
index=0;
[T R1]=ode45(@Repelling1,timeR,R0,opts);
clear Eigencheck;
index=0;
[T R2]=ode45(@Repelling2,timeR,R0,opts);
clear Eigencheck;
index=0;
%scatter(Y(:,1),Y(:,2))

figure
hold on
quiver(x,y,u,v,'k')
plot(A1(:,1),A1(:,2),'b','LineWidth',2)
plot(A2(:,1),A2(:,2),'b','LineWidth',2)
plot(R1(:,1),R1(:,2),'r','LineWidth',2)
plot(R2(:,1),R2(:,2),'r','LineWidth',2)
scatter(s2xmax, s2ymax,'m','filled')
scatter(s1xmin(1), s1ymin(1),'m','filled')
%contour(x,y,abs(Eigenvalue1))
%figure
%quiver(x,y,Eigenvector1(:,:,1),Eigenvector1(:,:,2),'k')
%}



function dzdt=Attracting1(t,z)
    index=index+1;
    Eigencheck(index)=abs(interp2(x,y,Eigenvalue1,z(1),z(2)));
    dzdt = zeros(3,1);
    if and(and(z(1)<=xmax,z(1)>=xmin),and(z(2)<=ymax,z(2)>=ymin))
        [Ix2 Iy2]=smooth_xi(z(1),z(2),x,y,Eigenvector2);
        dzdt(1) = Ix2; 
        dzdt(2) = Iy2;
        if index > 1
            if Eigencheck(index-1)>=Eigencheck(index)
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
end

function dzdt=Attracting2(t,z)
    index=index+1;
    Eigencheck(index)=abs(interp2(x,y,Eigenvalue1,z(1),z(2)));
    dzdt = zeros(3,1);
    if and(and(z(1)<=xmax,z(1)>=xmin),and(z(2)<=ymax,z(2)>=ymin))
        [Ix2 Iy2]=smooth_xi(z(1),z(2),x,y,Eigenvector2);
        dzdt(1) = -Ix2; 
        dzdt(2) = -Iy2;
        if index > 1
            if Eigencheck(index-1)>=Eigencheck(index)
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
end

function dzdt=Repelling1(t,z)
    index=index+1;
    Eigencheck(index)=abs(interp2(x,y,Eigenvalue2,z(1),z(2)));
    dzdt = zeros(3,1);
    if and(and(z(1)<=xmax,z(1)>=xmin),and(z(2)<=ymax,z(2)>=ymin))
        [Ix1 Iy1]=smooth_xi(z(1),z(2),x,y,Eigenvector1);
        dzdt(1) = Ix1; 
        dzdt(2) = Iy1;
        if index > 1
            if Eigencheck(index-1)>=Eigencheck(index)
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
end

function dzdt=Repelling2(t,z)
    index=index+1;
    Eigencheck(index)=abs(interp2(x,y,Eigenvalue2,z(1),z(2)));
    dzdt = zeros(3,1);
    if and(and(z(1)<=xmax,z(1)>=xmin),and(z(2)<=ymax,z(2)>=ymin))
        [Ix1 Iy1]=smooth_xi(z(1),z(2),x,y,Eigenvector1);
        dzdt(1) = -Ix1; 
        dzdt(2) = -Iy1;
        if index > 1
            if Eigencheck(index-1)>=Eigencheck(index)
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
end








function [value,isterminal,direction]=eventfun(t,z)
    value=[z(3)-eps,z(3)+eps];
    isterminal=[1,1];
    direction=[0,0];
end
%}


%{
Ax0=s2xmax(1);
Ay0=s2ymax(1);
A0=[Ax0,Ay0,Ax0,Ay0];
Rx0=s1xmin(1);
Ry0=s1ymin(1);
R0=[Rx0,Ry0,Rx0,Ry0];
timeA=[0 0.5];
timeR=[0 0.95];
%EigcheckA
%EigcheckA
%opts=odeset('RelTol',1e-9,'Events', @eventfun);
opts=odeset('Events', @eventfun);
[T A]=ode45(@Attracting,timeA,A0,opts);
[T R]=ode45(@Repelling,timeR,R0,opts);
%scatter(Y(:,1),Y(:,2))
figure
hold on
quiver(x,y,u,v,'k')
plot(A(:,1),A(:,2),'b','LineWidth',2)
plot(A(:,3),A(:,4),'b','LineWidth',2)
plot(R(:,1),R(:,2),'r','LineWidth',2)
plot(R(:,3),R(:,4),'r','LineWidth',2)
scatter(s2xmax, s2ymax,'m','filled')
scatter(s1xmin(1), s1ymin(1),'m','filled')


function dzdt=Attracting(t,z)
    dzdt = zeros(4,1);
    [Ix1 Iy1]=smooth_xi(z(1),z(2),x,y,Eigenvector2);
    [Ix2 Iy2]=smooth_xi(z(3),z(4),x,y,Eigenvector2);
    dzdt(1) = Ix1; 
    dzdt(2) = Iy1; 
    dzdt(3) = -Ix2; 
    dzdt(4) = -Iy2; 
        
end


function dzdt=Repelling(t,z)
    dzdt = zeros(4,1);
    [Ix1 Iy1]=smooth_xi(z(1),z(2),x,y,Eigenvector1);
    [Ix2 Iy2]=smooth_xi(z(3),z(4),x,y,Eigenvector1);
    dzdt(1) = Ix1; 
    dzdt(2) = Iy1; 
    dzdt(3) = -Ix2; 
    dzdt(4) = -Iy2; 
        
end


function [value,isterminal,direction]=eventfun(t,y)
    value=1;
    isterminal=1;
    direction=0;
end


%}


%{
for i=1:1:ylength
    for j=1:1:xlength
        for k=1:1:length(s1xmin)
            if and(abs(s1xmin(k)-x(i,j))<1e-5,abs(s1ymin(k)-y(i,j))<1e-5)
               quiver(x(i,j),y(i,j),1/2*Eigenvector(i,j,2,1),1/2*Eigenvector(i,j,2,2),'r','linewidth',2)
               quiver(x(i,j),y(i,j),-1/2*Eigenvector(i,j,2,1),-1/2*Eigenvector(i,j,2,2),'r','linewidth',2)
            end
        end
    end
end
%}    



%T=eta_tracing(s1xmin,s1ymin,linspace(0,0.1,3),x1,y1,Eigenvalue(:,:,1),Eigenvalue(:,:,2),Eigenvector(:,:,1,:),Eigenvector(:,:,2,:),1,1,odeset('Reltol',1e-5))
%[p t]=eta_tracing(s1xmin(1),s1ymin(1),linspace(0,0.1,10),x1,y1,Eigenvalue(:,:,1),Eigenvalue(:,:,2),Eigenvector(:,:,1,:),Eigenvector(:,:,2,:),-1,0.1,odeset('Reltol',1e-5))
%plot(p,t,'k*')


%[T Y]=ode113(@odefun,[0 10], [1 1])
%plot(Y(:,1), Y(:,2) )


%function dydt = odefun(t,y)
%    dydt = zeros(2,1);
%    dydt(1) = y(1)
%    dydt(2) = y(2)
    

%{

function dzdt=Attracting2(t,z)
    index=index+1;
    Eigencheck(index)=abs(interp2(x,y,Eigenvalue1,z(1),z(2)));
    dzdt = zeros(3,1);
    [Ix2 Iy2]=smooth_xi(z(1),z(2),x,y,Eigenvector2);
    dzdt(1) = -Ix2; 
    dzdt(2) = -Iy2;
    if index > 1
        if Eigencheck(index-1)>=Eigencheck(index)
            dzdt(3)=0;
        else
            dzdt(3)=1;
        end
    end
end

function dzdt=Repelling1(t,z)
    index=index+1;
    Eigencheck(index)=abs(interp2(x,y,Eigenvalue2,z(1),z(2)));
    dzdt = zeros(3,1);
    [Ix1 Iy1]=smooth_xi(z(1),z(2),x,y,Eigenvector1);
    dzdt(1) = Ix1; 
    dzdt(2) = Iy1; 
    if index > 1
        if Eigencheck(index-1)>=Eigencheck(index)
            dzdt(3)=0;
        else
            dzdt(3)=1;
        end
    end
end

function dzdt=Repelling2(t,z)
    index=index+1;
    Eigencheck(index)=abs(interp2(x,y,Eigenvalue2,z(1),z(2)));
    dzdt = zeros(3,1);
    [Ix1 Iy1]=smooth_xi(z(1),z(2),x,y,Eigenvector1);
    dzdt(1) = -Ix1;
    dzdt(2) = -Iy1; 
    if index > 1
        if Eigencheck(index-1)>=Eigencheck(index)
            dzdt(3)=0;
        else
            dzdt(3)=1;
        end
    end
end
%}

end