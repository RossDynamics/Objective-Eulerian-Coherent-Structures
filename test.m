function test()

%Input Parameters
xmin=-0;
xmax=2;
ymin=-0.5;
ymax=0.5;
res=10;
dist=1/res;
xlength=((xmax-xmin)*res+1);
ylength=((ymax-ymin)*res+1);

global x y u v Eigenvalue1 Eigenvalue2 Eigenvector1 Eigenvector2

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

%quiver(x,y,Eigenvector1(:,:,1),Eigenvector1(:,:,2))
[x1 y1]=smooth_xi(1.25,0,x,y,Eigenvector1)
[x2 y2]=smooth_xi(1.25,0,x,y,Eigenvector2)



end