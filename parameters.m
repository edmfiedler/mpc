%% Parameters

% Outplet pipe cm^2
a1 = 1.2272;
a2 = 1.2272;
a3 = 1.2272;
a4 = 1.2272;
a = [a1;a2;a3;a4];

% Cross-sectional area cm^2
A1 = 380.1327;
A2 = 380.1327;
A3 = 380.1327;
A4 = 380.1327;
A = [A1;A2;A3;A4];

% Flow distribution
gam1 = 0.6;
gam2 = 0.75;
gam = [gam1;gam2];

% Acceleration of gravity cms^-2
g = 981;

% Density of water gcm^-3
rho = 1;

p = [a;A;gam;g;rho];