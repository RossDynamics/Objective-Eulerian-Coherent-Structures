function [Eigenvalue1 Eigenvalue2 Eigenvector1 Eigenvector2] = CalculateEigenFields(S)
%CALCULATEEIGENFIELDS Calculates Eigenvalue and Eigenvector Fields
%   Calculates Eigenvalue and Eigenvector Fields from the Eulerian Strain
%   Tensor, S, for a 2D flow.

xlength=length(S(1,:,1,1));
ylength=length(S(:,1,1,1));

for i=[1:1:ylength]
    for j=[1:1:xlength]
        Matrix(1,1)=S(i,j,1,1);
        Matrix(1,2)=S(i,j,1,2);
        Matrix(2,1)=S(i,j,2,1);
        Matrix(2,2)=S(i,j,2,2);
        Eigenvalues=eig(Matrix);
        [Eigenvectors D]=eig(Matrix);
        if Eigenvalues(1)<Eigenvalues(2)
            Eigenvalue1(i,j)=Eigenvalues(1);
            Eigenvalue2(i,j)=Eigenvalues(2);
            Eigenvector1(i,j,:)=Eigenvectors(:,1);
            Eigenvector2(i,j,:)=Eigenvectors(:,2);
        else
            Eigenvalue1(i,j)=Eigenvalues(2);
            Eigenvalue2(i,j)=Eigenvalues(1);
            Eigenvector1(i,j,:)=Eigenvectors(:,2);
            Eigenvector2(i,j,:)=Eigenvectors(:,1);
        end
        clear Matrix Eigenvalues Eigenvectors D
    end
end


end

