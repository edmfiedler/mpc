function [xs,hs] = findsteadystate(u,d,p)

xs0 = 5000*ones(4,1);
us = u;
ds = d;
xs = fsolve(@ModifiedFourTankSystemWrap,xs0,[],us,ds,p);
hs = xs./(p(12).*p(5:8));

end

