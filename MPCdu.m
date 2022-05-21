function DU = MPCdu(r,n,u)

DU = zeros(n*size(r,1),1);
for i = 1:size(r,1):size(DU,1)
    DU(i:i+size(r,1)-1,1) = [r(1);r(2)];
end

DU(1:size(r,1)) = DU(1:size(r,1))+u;

end

