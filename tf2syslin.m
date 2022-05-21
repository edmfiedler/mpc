function sysd = tf2syslin(TF,Ts);
addpath('Realization')

si = length(TF)/2;

% Allocate memory
num = cell(si,si);
den = cell(si,si);
lambda = zeros(si,si);

% Build Array
[n1,d1] = tfdata(TF(1));
[n2,d2] = tfdata(TF(2));
[n3,d3] = tfdata(TF(3));
[n4,d4] = tfdata(TF(4));

num(1,1) = n1;
den(1,1) = d1;

num(1,2) = n2;
den(1,2) = d2;

num(2,1) = n3;
den(2,1) = d3;

num(2,2) = n4;
den(2,2) = d4;

Nmax = 2*1200;
tol = 1e-8;

[Ad,Bd,Cd,Dd,sH] = mimoctf2dss(num,den,lambda,Ts,Nmax,tol);

sysd = ss(Ad,Bd,Cd,Dd,Ts);
end

