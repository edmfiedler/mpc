function C = eco_cost(seed,len)
randn('state',seed)

C = zeros(len,1);
c1 = randn(1,length(C)/2);
c2 = c1;
C(1:2:end) = c1;
C(2:2:end) = c2;
C = (C - min(C));
for i = 1:2:length(C)
   if i < 150
        C(i) = C(i)+10;
   end
end

for i = 2:2:length(C)
   if i > 250
        C(i) = C(i)+10;
   end
end

end

