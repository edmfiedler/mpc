clear
clc

%% Load simulation parameters
parameters

%% Determine steady state from which to step from

% input selection cm^3s^-1 
F1 = 250;
F2 = 325;
F3 = 100;
F4 = 100;

u = [F1;F2];
d = [F3;F4];

[xs,hs] = findsteadystate(u,d,p);

%% Apply linmod to simulink model to obtain state-space representation

sys_l = linmod('systemlin',xs,[u;d]);

A = sys_l.a;
B = sys_l.b(:,1:2);
Bv = sys_l.b(:,3:4);
C = sys_l.c;
D = zeros(2,2);

%% Obtain the transfer function representation

[n1,d1] = ss2tf(A,B,C,D,1);
[n2,d2] = ss2tf(A,B,C,D,2);

tf_1 = minreal(tf(n1(1,:),d1));
tf_2 = minreal(tf(n2(1,:),d2));
tf_3 = minreal(tf(n1(2,:),d1));
tf_4 = minreal(tf(n2(2,:),d2));

%% Obtain the discrete-time state-space description

% Same Ts as used in deterministic_step, 9s
Ts = 5;

sysc = ss(A,B,C,D);
sysd = c2d(sysc,Ts);

Ad = sysd.a;
Bd = sysd.b;
Cd = sysd.c;
Dd = sysd.d;
[~,Bvd] = c2d(A,Bv,Ts);

%% Obtain Markov Parameters

n = floor(Ts)*100;

H_lin = markovvec(Ad,Bd,Cd,n);

%% Compare with system from approximated transfer functiosn
xmax = 400;

load('H_approx.mat') % Obtained from deterministic_step

figure('Units','inches','Position',[0 0 6.693 4],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

subplot(2,2,1)
plot(0:n-1,H_lin(1:2:end,1),'k','linewidth',1.1)
hold on
plot(0:n-1,H(1:2:end,1),'k:','linewidth',1.1)
hold off
grid('on')
xlim([0 xmax])
legend('Linearised','Approximated t.f.')
%title('$u_1$ $\to$ $y_1$','interpreter','latex','FontSize',12)
xlabel('Step [k]')

subplot(2,2,2)
plot(0:n-1,H_lin(1:2:end,2),'k','linewidth',1.1)
hold on
plot(0:n-1,H(1:2:end,2),'k:','linewidth',1.1)
hold off
grid('on')
xlim([0 xmax])
legend('Linearised','Approximated t.f.')
%title('$u_2$ $\to$ $y_1$','interpreter','latex','FontSize',12)
xlabel('Step [k]')

subplot(2,2,3)
plot(0:n-1,H_lin(2:2:end,1),'k','linewidth',1.1)
hold on
plot(0:n-1,H(2:2:end,1),'k:','linewidth',1.1)
hold off
grid('on')
xlim([0 xmax])
legend('Linearised','Approximated t.f.')
%title('$u_1$ $\to$ $y_2$','interpreter','latex','FontSize',12)
xlabel('Step [k]')

subplot(2,2,4)
plot(0:n-1,H_lin(2:2:end,2),'k','linewidth',1.1)
hold on
plot(0:n-1,H(2:2:end,2),'k:','linewidth',1.1)
hold off
grid('on')
xlim([0 xmax])
legend('Linearised','Approximated t.f.')
%title('$u_2$ $\to$ $y_2$','interpreter','latex','FontSize',12)
xlabel('Step [k]')

sgtitle('Comparison of Markov Parameters','FontSize',11)

set(subplot(2,2,1), 'Position', [.055, .55, .39, .35]);
set(subplot(2,2,2), 'Position', [.555, .55, .39, .35]);
set(subplot(2,2,3), 'Position', [.055, .085, .39, .35]);
set(subplot(2,2,4), 'Position', [.555, .085, .39, .35]);




