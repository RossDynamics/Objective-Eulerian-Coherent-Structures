function [type] = ClassifySingularity(eigenvector)
%CLASSIFYSINGULARITIES Summary of this function goes here
%   Detailed explanation goes here
% 3 if trisector, 2 if wedge, 0 otherwise

%Matrix = [-1 1;1 +1];
%[V D] = eig(Matrix);
%eigenvector=V(:,1);

tol=1e-5;
theta=[0:0.001:2*pi];
radialvector=[cos(theta);sin(theta)];
 for i=1:1:length(theta)
     S=singularityProduct(eigenvector, radialvector(:,i));
     if S>1-tol
         F(i)=1;
     elseif S<tol
         F(i)=0;
     else
        F(i)=2;
     end
 end

 zero=0;
 one=0;
 if F(1)==0
     zero=zero+1;
 elseif F(1)==1
     one=one+1;
 end
 
 for i=2:1:length(F)
    if ne(F(i),F(i-1))
        if F(i)==0
            zero=zero+1;   
        elseif F(i)==1
            one=one+1;
        end
    end
 end

 if and(one==3,zero==3)
     type = 3;
 elseif and(one==3,zero==1)
     type = 2;
 else
     type = 0;
 end
 
 
end

function f=singularityProduct(eigenvect, radialVector)
    f = abs(dot(eigenvect,radialVector))/(norm(eigenvect)*norm(radialVector));
end