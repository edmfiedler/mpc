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

zd = []; Td = []; ud = []; ref = []; x = xs; xd = []; u = [F1;F2]; d = [F3;F4]; ds = d;
xh = [0;0;0;0]; Pk = 10000000; xhd =[]; Pd = []; yd = [];
randn('state',1200)
for k = k_0:k_f
    %Pk = Pss; % Static
    
    % Simulation
    if k < 500
        d = [F3;F4]+chol(Qk)'*randn(2,1);
    else 
        d = [F3;F4]+chol(Qk)'*randn(2,1)+[25;25];
    end
    
    % Simulation
    [T,X] = ode15s(@ModifiedFourTankSystem,[k*Ts (k+1)*Ts],x,[],u,d,p);
    x = real(X(end,:)');
    xd = [xd x];
    Td = [Td;T(1:end-1)];
    y = Cd*x+chol(Rk)'*randn(2,1);
    yd = [yd;y'];
    ud = [ud; u' d'];
    
    yk = y-hs(1:2);
    uk = u-[F1;F2];
    % State Estimation
    Pk = Pss;
    [xh,Pk] = statestimator(xh,Pk,Ad,Bd,Bvd,Cd,uk,yk,Rk,Qk);
    xhd = [xhd xh+xs];
    
    if k == 200
        u = [F1;F2]+[0;u_step];
    end
    
    if k == 600
        u = [F1;F2]+[u_step;u_step];
    end
    
    if k == 800
        u = [F1;F2]+[u_step;0];
    end
end

%% Plotting the normalized steps

figure('Units','inches','Position',[0 0 6.693 4],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

xl = 5000;

subplot(2,2,1)
stairs(k_0:Ts:k_f*Ts,xd(1,:),'k','linewidth',1)
hold on
stairs(k_0:Ts:k_f*Ts,xhd(1,:),'k--','linewidth',1)
hold off
ax = gca;
ax.YAxis.Exponent = 4;
grid('on')
xlim([0 xl])
ylim([1.3e4 1.8e4])
ylabel('$x_1$ [g]')
xlabel('Time [s]')
legend('$x$','$\hat{x}$','Location','NorthWest')


subplot(2,2,2)
stairs(k_0:Ts:k_f*Ts,xd(2,:),'k','linewidth',1)
hold on
stairs(k_0:Ts:k_f*Ts,xhd(2,:),'k--','linewidth',1)
hold off
ax = gca;
ax.YAxis.Exponent = 4;
grid('on')
xlim([0 xl])
ylim([2.4e4 3.05e4])
ylabel('$x_2$ [g]')
xlabel('Time [s]')
legend('$x$','$\hat{x}$','Location','NorthWest')

subplot(2,2,3)
stairs(k_0:Ts:k_f*Ts,xd(3,:),'k','linewidth',1)
hold on
stairs(k_0:Ts:k_f*Ts,xhd(3,:),'k--','linewidth',1)
hold off
ax = gca;
ax.YAxis.Exponent = 4;
grid('on')
xlim([0 xl])
ylim([0.38e4 0.6e4])
ylabel('$x_3$ [g]')
xlabel('Time [s]')
legend('$x$','$\hat{x}$','Location','NorthWest')

subplot(2,2,4)
stairs(k_0:Ts:k_f*Ts,xd(4,:),'k','linewidth',1)
hold on
stairs(k_0:Ts:k_f*Ts,xhd(4,:),'k--','linewidth',1)
hold off
ax = gca;
ax.YAxis.Exponent = 4;
grid('on')
xlim([0 xl])
ylim([0.45e4 0.74e4])
ylabel('$x_4$ [g]')
xlabel('Time [s]')
legend('$x$','$\hat{x}$','Location','NorthWest')

sgtitle('Static Kalman filtering with process noise step change at 2500s','FontSize',11)

set(subplot(2,2,1), 'Position', [.055, .55, .39, .35]);
set(subplot(2,2,2), 'Position', [.555, .55, .39, .35]);
set(subplot(2,2,3), 'Position', [.055, .085, .39, .35]);
set(subplot(2,2,4), 'Position', [.555, .085, .39, .35]);