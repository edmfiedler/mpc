function H = markovvec(A,B,C,n)

si = size(C,1);
H = zeros(si*n,si);

% Intermediate variable
H_i = C;
for i = 1:si:n*si
    H(i:i+si-1,:) = H_i*B;
    H_i = H_i*A;
end

end

