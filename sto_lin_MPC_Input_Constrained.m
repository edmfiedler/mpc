clear
clc

%% Obtain the discrete-time state-space system

parameters

F1 = 250;
F2 = 325;
F3 = 100;
F4 = 100;

u = [F1;F2];
d = [F3;F4];

[xs,hs] = findsteadystate(u,d,p);

sys_l = linmod('systemlin',xs,[u;d]);

A = sys_l.a;
B = sys_l.b(:,1:2);
Bv = sys_l.b(:,3:4);
C = sys_l.c;
D = zeros(2,2);

% Sample time
Ts = 5;

sysc = ss(A,B,C,D);
sysd = c2d(sysc,Ts);

Ad = sysd.a;
Bd = sysd.b;
Cd = sysd.c;
Dd = sysd.d;
[~,Bvd] = c2d(A,Bv,Ts);

%% Compute the MPC related matrices

% Depth
n = 30*Ts;

% Weighting Matrices
Q = [1 0;0 1];
S = [1 0;0 1];

% Build Matrices for Control Law
[H,Mx0,Mr,Mu_1,Md,Aqp,Phi,Gamma_d] = MPCmat(Ad,Bd,Cd,Bvd,Q,S,n);
% For input constraint, only need Lambda from Aqp
Aqp = Aqp(1:length(H),:);

%% Simulate with MPC

% Simulation length
k_0 = 0;
k_f = 250;
sim_length = k_0:k_f;

% Create a ramp and step in the simulation
ch_ref = ceil(k_f/3);
r1 = [(hs(1):((85-hs(1))/ch_ref):85) 85*ones(1,ch_ref) 10*ones(1,ch_ref+n)]'-hs(1);
r2 = [(hs(2):((85-hs(2))/ch_ref):85) 85*ones(1,ch_ref) 10*ones(1,ch_ref+n)]'-hs(2);
ref = zeros(length(r1)+length(r2),1);
ref(1:2:end) = r1;
ref(2:2:end) = r2;

% Set upper and lower limits for Inputs
u_min = [0;0]-[F1;F2];
u_max = [600;600]-[F1;F2];
U_min = repmat(u_min,n,1);
U_max = repmat(u_max,n,1);

% Set random variables
randn('state',1200)
r_and = randn(2,2*length(sim_length));
d_rand = r_and(:,1:0.5*length(r_and));
v_rand = r_and(:,0.5*length(r_and)+1:end);

% Initial conditions
zd = []; yd = []; Td = []; ud = []; set_p = []; x = [0;0;0;0]; u = [0;0]; d = [0;0];
xh = [0;0;0;0]; Pk = 10; xhd =[]; Pd = []; Qk = 100*eye(2); Rk = 1*eye(2);
for k = sim_length
    % Simulation
    d = chol(Qk)'*d_rand(:,k+1);
    x = Ad*x+Bd*u+Bvd*d;
    v = chol(Rk)'*v_rand(:,k+1);
    y = Cd*x+v;
    yd = [yd;y'];
    z = Cd*x;
    zd = [zd;z'];
    ud = [ud; u' d'];
    
    % State Estimation
    [xh,Pk] = statestimator(xh,Pk,Ad,Bd,Bvd,Cd,u,y,Rk,Qk);
    xhd = [xhd xh];
    
    % Select reference
    i_ref = 2*k+1;
    R = ref(i_ref:(i_ref+2*n-1));
    % Control Law
        % Delta U
    d_u_min = [-10;-10];
    D_U_min = MPCdu(d_u_min,n,u);
    d_u_max = [10;10];
    D_U_max = MPCdu(d_u_max,n,u);
        % solving QP
    gu = Mx0*xh+Mr*R+Mu_1*u;
    [u,info] = qpsolver(H,gu,U_min,U_max,Aqp,[D_U_min],[D_U_max],[]);
    u = [u(1);u(2)];
    
    set_p = [set_p;R(1:2)'+hs(1:2)'];
end

%% Plot outputs and inputs

figure('Units','inches','Position',[0 0 6.693 4],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

subplot(2,1,1)
stairs(k_0:Ts:k_f*Ts,yd(:,1)+hs(1),'k','linewidth',1.3)
hold on
stairs(k_0:Ts:k_f*Ts,yd(:,2)+hs(2),'Color', [0.5 0.5 0.5],'linewidth',1.3)
stairs(k_0:Ts:k_f*Ts,set_p(:,1),'k:','linewidth',0.8)
stairs(k_0:Ts:k_f*Ts,set_p(:,2),':','Color', [0.5 0.5 0.5],'linewidth',0.8)
hold off
grid('on')
xlabel('Time [s]')
ylabel('$y$ [cm]')
xlim([k_0 k_f*Ts])
ylim([0 100])
legend('$y_1$','$y_2$','$r_1$','$r_2$')

subplot(2,1,2)
stairs(k_0:Ts:k_f*Ts,ud(:,1)+F1,'k','linewidth',1.3)
hold on
stairs(k_0:Ts:k_f*Ts,ud(:,2)+F2,'Color', [0.5 0.5 0.5],'linewidth',1.3)
hold off
grid('on')
xlabel('Time [s]')
ylabel('$u$ [cm$^3$ s$^{-1}$]')
xlim([k_0 k_f*Ts])
ylim([-80 700])
legend('$u_1$','$u_2$')

sgtitle('Input-constrained MPC controlling linearised stochastic model','FontSize',11)

set(subplot(2,1,1), 'Position', [.065, .55, .90, .35]);
set(subplot(2,1,2), 'Position', [.065, .085, .90, .35]);