[a,b]=size(x);
index=1
q=zeros(size(x));
w=zeros(size(x));

for i=[1:1:a]
    for j=[1:1:b]
       if mod(index,4)==1
           q(i,j)=x(i,j);
           w(i,j)=y(i,j);
       else
           q(i,j)=NaN;
           w(i,j)=NaN;
       end
       index=index+1;
    end
    index=index+2;
end
clear a b i j index