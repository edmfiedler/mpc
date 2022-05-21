function R = MPCref(r,n)

R = zeros(n*size(r,1),1);
for i = 1:size(r,1):size(R,1)
    R(i:i+size(r,1)-1,1) = [r(1);r(2)];
end

end

