olon = X(1,1);
olat = Y(1,1);
X=dx
Y=dy
x=zeros(size(X));
y=zeros(size(X));
for i=[1:1:length(X(:,1))]
    for j=[1:1:length(X(1,:))]
        slon = X(i,j);
        slat = Y(i,j);
        [kmx kmy]=ll2kmxy(slat,slon,olat,olon);
        x(i,j)=kmx;
        y(i,j)=kmy;
    end 
end
clear i j kmx kmy olat olon slat slon