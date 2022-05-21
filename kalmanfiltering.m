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

%load('sys_approx.mat')

Ad = sysd.a;
Bd = sysd.b;
Cd = sysd.c;
Dd = sysd.d;
[~,Bvd] = c2d(A,Bv,Ts);

%% Simulate

% Simulation length
k_0 = 0;
k_f = 1000;

u_step = 10;

Qk = 100*eye(2); Rk = 10*eye(2);
[~,~,Pss] = kalman(ss(Ad,[Bd Bvd],Cd,zeros(2,4),Ts),Qk,Rk);

zd = []; Td = []; ud = []; ref = []; x = [0;0;0;0]; xd = []; u = [0;0]; d = [0;0]; ds = d;
xh = [50000;50000;50000;50000]; Pk = 10; xhd =[]; Pd = [];
randn('state',1200)
for k = k_0:k_f
    %Pk = Pss; % Static
    
    % Simulation
    if k < 500
        d = chol(Qk)'*randn(2,1);
    else 
        d = chol(Qk)'*randn(2,1)+[25;25];
    end
    %d = [0;0];
    x = Ad*x+Bd*u+Bvd*d;
    xd = [xd x];
    v = chol(Rk)'*randn(2,1);
    %v = [0;0];
    y = Cd*x+v;
    z = Cd*x;
    zd = [zd;z'];
    ud = [ud; u' d'];
    
    % State Estimation
    Pk = Pss;
    [xh,Pk] = statestimator(xh,Pk,Ad,Bd,Bvd,Cd,u,y,Rk,Qk);
    xhd = [xhd xh];
    
    if k == 200
        u = [0;u_step];
    end
    
    if k == 600
        u = [u_step;u_step];
    end
    
    if k == 800
        u = [u_step;0];
    end
end

%% Plotting the normalized steps

figure('Units','inches','Position',[0 0 6.693 4],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

xl = 5000;

subplot(2,2,1)
stairs(k_0:Ts:k_f*Ts,xd(1,:)+xs(1),'k','linewidth',1)
hold on
stairs(k_0:Ts:k_f*Ts,xhd(1,:)+xs(1),'k--','linewidth',1)
hold off
ax = gca;
ax.YAxis.Exponent = 4;
grid('on')
xlim([0 xl])
ylim([1e4 2e4])
ylabel('$x_1$ [g]')
xlabel('Time [s]')
legend('$x$','$\hat{x}$','Location','SouthEast')


subplot(2,2,2)
stairs(k_0:Ts:k_f*Ts,xd(2,:)+xs(2),'k','linewidth',1)
hold on
stairs(k_0:Ts:k_f*Ts,xhd(2,:)+xs(2),'k--','linewidth',1)
hold off
ax = gca;
ax.YAxis.Exponent = 4;
grid('on')
xlim([0 xl])
ylim([2e4 3e4])
ylabel('$x_2$ [g]')
xlabel('Time [s]')
legend('$x$','$\hat{x}$','Location','SouthEast')

subplot(2,2,3)
stairs(k_0:Ts:k_f*Ts,xd(3,:)+xs(3),'k','linewidth',1)
hold on
stairs(k_0:Ts:k_f*Ts,xhd(3,:)+xs(3),'k--','linewidth',1)
hold off
ax = gca;
ax.YAxis.Exponent = 4;
grid('on')
xlim([0 xl])
ylim([0.2e4 1.5e4])
ylabel('$x_3$ [g]')
xlabel('Time [s]')
legend('$x$','$\hat{x}$','Location','NorthEast')

subplot(2,2,4)
stairs(k_0:Ts:k_f*Ts,xd(4,:)+xs(4),'k','linewidth',1)
hold on
stairs(k_0:Ts:k_f*Ts,xhd(4,:)+xs(4),'k--','linewidth',1)
hold off
ax = gca;
ax.YAxis.Exponent = 4;
grid('on')
xlim([0 xl])
ylim([0.2e4 1.5e4])
ylabel('$x_4$ [g]')
xlabel('Time [s]')
legend('$x$','$\hat{x}$','Location','NorthEast')

sgtitle('Dynamic Kalman filtering with process noise step change at 2500s','FontSize',11)

set(subplot(2,2,1), 'Position', [.055, .55, .39, .35]);
set(subplot(2,2,2), 'Position', [.555, .55, .39, .35]);
set(subplot(2,2,3), 'Position', [.055, .085, .39, .35]);
set(subplot(2,2,4), 'Position', [.555, .085, .39, .35]);