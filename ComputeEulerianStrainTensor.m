function S = ComputeEulerianStrainTensor(dxu,dyu,dxv,dyv)
%COMPUTEEULERIANSTRAINTENSOR computers the eulerian strain tensor
%   computers the eulerian strain tensor from the gradient vectors of the
%   velocity fields
%   dxu = du/dx, dyu = du/dy
%   dxv = dv/dx, dyv = dv/dy
xlength=length(dxu(1,:));
ylength=length(dyv(:,1));

for i=[1:1:ylength]
    for j=[1:1:xlength]
        S(i,j,1,1)=1/2*(dxu(i,j)+dxu(i,j));
        S(i,j,1,2)=1/2*(dxv(i,j)+dyu(i,j));
        S(i,j,2,1)=1/2*(dyu(i,j)+dxv(i,j));
        S(i,j,2,2)=1/2*(dyv(i,j)+dyv(i,j));
    end
end

end

