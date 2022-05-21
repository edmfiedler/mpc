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

%% Simulation Scenario - Deterministic Step Response

% Simulation time frame
t_0 = 0;
t_f = 1300;
t_step = 0;

% Simulations
pch = 1.1;
[y_1_10,u_1_10,t_1_10] = detsim(pch,1,u,d,p,xs,t_0,t_f,t_step);
[tf_1_10_1,tf_1_10_2,K_1_10,tau_1_10] = firstordertf(y_1_10,u_1_10,t_1_10,hs,t_step);

[y_2_10,u_2_10,t_2_10] = detsim(pch,2,u,d,p,xs,t_0,t_f,t_step);
[tf_2_10_1,tf_2_10_2,K_2_10,tau_2_10] = firstordertf(y_2_10,u_2_10,t_2_10,hs,t_step);

pch = 1.25;
[y_1_25,u_1_25,t_1_25] = detsim(pch,1,u,d,p,xs,t_0,t_f,t_step);
[tf_1_25_1,tf_1_25_2,K_1_25,tau_1_25] = firstordertf(y_1_25,u_1_25,t_1_25,hs,t_step);

[y_2_25,u_2_25,t_2_25] = detsim(pch,2,u,d,p,xs,t_0,t_f,t_step);
[tf_2_25_1,tf_2_25_2,K_2_25,tau_2_25] = firstordertf(y_2_25,u_2_25,t_2_25,hs,t_step);

pch = 1.5;
[y_1_50,u_1_50,t_1_50] = detsim(pch,1,u,d,p,xs,t_0,t_f,t_step);
[tf_1_50_1,tf_1_50_2,K_1_50,tau_1_50] = firstordertf(y_1_50,u_1_50,t_1_50,hs,t_step);

[y_2_50,u_2_50,t_2_50] = detsim(pch,2,u,d,p,xs,t_0,t_f,t_step);
[tf_2_50_1,tf_2_50_2,K_2_50,tau_2_50] = firstordertf(y_2_50,u_2_50,t_2_50,hs,t_step);
%% Plotting steps

close all
figure('Units','inches','Position',[0 0 6.693 4],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');


subplot(2,2,1)
plot(t_1_10,y_1_10(:,1),'k','linewidth',1.1)
hold on
plot(t_1_25,y_1_25(:,1),'Color', [0.4 0.4 0.4],'linewidth',1.1)
plot(t_1_50,y_1_50(:,1),'Color', [0.7 0.7 0.7],'linewidth',1.1)
hold off
xlim([t_0 t_f])
grid('on')
ylabel('$y_1$ [cm]')
xlabel('Time [s]')
ylim([35 60])
title('Step change in $u_1$','FontSize',10)
legend('10\%','25\%','50\%','Location','East')

subplot(2,2,3)
plot(t_1_10,y_1_10(:,2),'k','linewidth',1.1)
hold on
plot(t_1_25,y_1_25(:,2),'Color', [0.4 0.4 0.4],'linewidth',1.1)
plot(t_1_50,y_1_50(:,2),'Color', [0.7 0.7 0.7],'linewidth',1.1)
hold off
xlim([t_0 t_f])
ylim([65 85])
grid('on')
ylabel('$y_2$ [cm]')
xlabel('Time [s]')
legend('10\%','25\%','50\%','Location','East')

subplot(2,2,2)
plot(t_2_10,y_2_10(:,1),'k','linewidth',1.1)
hold on
plot(t_2_25,y_2_25(:,1),'Color', [0.4 0.4 0.4],'linewidth',1.1)
plot(t_2_50,y_2_50(:,1),'Color', [0.7 0.7 0.7],'linewidth',1.1)
hold off
xlim([t_0 t_f])
ylim([35 50])
grid('on')
ylabel('$y_1$ [cm]')
xlabel('Time [s]')
title('Step change in $u_2$','FontSize',10)
legend('10\%','25\%','50\%','Location','East')

subplot(2,2,4)
plot(t_2_10,y_2_10(:,2),'k','linewidth',1.1)
hold on
plot(t_2_25,y_2_25(:,2),'Color', [0.4 0.4 0.4],'linewidth',1.1)
plot(t_2_50,y_2_50(:,2),'Color', [0.7 0.7 0.7],'linewidth',1.1)
hold off
xlim([t_0 t_f])
ylim([65 110])
grid('on')
ylabel('$y_{2}$ [cm]')
xlabel('Time [s]')
legend('10\%','25\%','50\%','Location','East')

sgtitle('Deterministic response to percentage changes in input from steady state','FontSize',11)

set(subplot(2,2,1), 'Position', [.055, .53, .39, .35]);
set(subplot(2,2,2), 'Position', [.555, .53, .39, .35]);
set(subplot(2,2,3), 'Position', [.055, .085, .39, .35]);
set(subplot(2,2,4), 'Position', [.555, .085, .39, .35]);

%% Plotting the normalized steps

figure('Units','inches','Position',[0 0 6.693 4],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');


subplot(2,2,1)
plot(t_1_10,y_1_10(:,1)./hs(1),'k','linewidth',1.1)
hold on
plot(t_1_25,y_1_25(:,1)./hs(1),'Color', [0.4 0.4 0.4],'linewidth',1.1)
plot(t_1_50,y_1_50(:,1)./hs(1),'Color', [0.7 0.7 0.7],'linewidth',1.1)
hold off
xlim([t_0 t_f])
grid('on')
ylabel('$y_1$ (normalized)')
xlabel('Time [s]')
title('Step change in $u_1$','FontSize',10)
legend('10\%','25\%','50\%','Location','East')

subplot(2,2,3)
plot(t_1_10,y_1_10(:,2)./hs(2),'k','linewidth',1.1)
hold on
plot(t_1_25,y_1_25(:,2)./hs(2),'Color', [0.4 0.4 0.4],'linewidth',1.1)
plot(t_1_50,y_1_50(:,2)./hs(2),'Color', [0.7 0.7 0.7],'linewidth',1.1)
hold off
xlim([t_0 t_f])
ylim([1 1.25])
grid('on')
ylabel('$y_2$ (normalized)')
xlabel('Time [s]')
legend('10\%','25\%','50\%','Location','East')

subplot(2,2,2)
plot(t_2_10,y_2_10(:,1)./hs(1),'k','linewidth',1.1)
hold on
plot(t_2_25,y_2_25(:,1)./hs(1),'Color', [0.4 0.4 0.4],'linewidth',1.1)
plot(t_2_50,y_2_50(:,1)./hs(1),'Color', [0.7 0.7 0.7],'linewidth',1.1)
hold off
xlim([t_0 t_f])
grid('on')
ylabel('$y_1$ (normalized)')
xlabel('Time [s]')
title('Step change in $u_2$','FontSize',10)
legend('10\%','25\%','50\%','Location','East')

subplot(2,2,4)
plot(t_2_10,y_2_10(:,2)./hs(2),'k','linewidth',1.1)
hold on
plot(t_2_25,y_2_25(:,2)./hs(2),'Color', [0.4 0.4 0.4],'linewidth',1.1)
plot(t_2_50,y_2_50(:,2)./hs(2),'Color', [0.7 0.7 0.7],'linewidth',1.1)
hold off
xlim([t_0 t_f])
ylim([1 1.7])
grid('on')
ylabel('$y_{2}$ (normalized)')
xlabel('Time [s]')
legend('10\%','25\%','50\%','Location','East')

sgtitle('Normalized deterministic response to percentage changes in input from steady state','FontSize',11)

set(subplot(2,2,1), 'Position', [.055, .53, .39, .35]);
set(subplot(2,2,2), 'Position', [.555, .53, .39, .35]);
set(subplot(2,2,3), 'Position', [.055, .085, .39, .35]);
set(subplot(2,2,4), 'Position', [.555, .085, .39, .35]);

%% Obtaining "Average" Transfer function

K_1_1 = mean([K_1_10(1) K_1_25(1) K_1_50(1)]);
K_1_2 = mean([K_1_10(2) K_1_25(2) K_1_50(2)]);

K_2_1 = mean([K_2_10(1) K_2_25(1) K_2_50(1)]);
K_2_2 = mean([K_2_10(2) K_2_25(2) K_2_50(2)]);

tau_1_1 = mean([tau_1_10(1) tau_1_25(1) tau_1_50(1)]);
tau_1_2 = mean([tau_1_10(2) tau_1_25(2) tau_1_50(2)]);

tau_2_1 = mean([tau_2_10(1) tau_2_25(1) tau_2_50(1)]);
tau_2_2 = mean([tau_2_10(2) tau_2_25(2) tau_2_50(2)]);

tf_1 = tf(K_1_1,[tau_1_1 1]);
tf_2 = tf(K_2_1,[tau_2_1 1]);
tf_3 = tf(K_1_2,[tau_1_2 1]);
tf_4 = tf(K_2_2,[tau_2_2 1]);

%% Obtain linear model from transfer function estimates

% Define sample time based on smallest tau
Ts = floor(min([tau_1_1 tau_1_2 tau_2_1 tau_2_2])/10);
Ts = 5;

sysd = tf2syslin([tf_1 tf_2 tf_3 tf_4],Ts);

%% Check validity by comparing to 25% step in u1 of noisy response
t = 0:Ts:1260;
u = [F1*1.25-F1;0].*ones(2,length(t));
[y,t] = lsim(sysd,u,t);

load('n_y_1_25.mat')
load('n_t_1_25.mat')

figure('Units','inches','Position',[0 0 6.693 2],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

subplot(1,2,1)
plot(t_1_25,y_1_25(:,1),':','Color', [0.4 0.4 0.4],'linewidth',1.1)
hold on
stairs(t,y(:,1)+hs(1),'k','linewidth',1.1)
hold off
xlim([0 1200])
grid('on')
ylabel('$y_1$ [cm]')
xlabel('Time [s]')
legend('Noisy response','Approximated response','Location','SouthEast')

subplot(1,2,2)
plot(t_1_25,y_1_25(:,2),':','Color', [0.4 0.4 0.4],'linewidth',1.1)
hold on
stairs(t,y(:,2)+hs(2),'k','linewidth',1.1)
hold off
xlim([0 1200])
grid('on')
ylabel('$y_2$ [cm]')
xlabel('Time [s]')
legend('Noisy response','Approximated response','Location','SouthEast')

sgtitle('Comparison of medium noisy response to approximated system','FontSize',11)

%% Compute Markov Parameters of approximated system

% Define depth
n = floor(Ts)*100;

H = markovvec(sysd.a,sysd.b,sysd.c,n);

%% Plot Markov Parameters

figure('Units','inches','Position',[0 0 6.693 4],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

n = 400;

subplot(2,2,1)
plot(H(1:2:end,1),'k','linewidth',1.1)
ax = gca;
ax.YAxis.Exponent = -3;
grid('on')
xlim([0 n])
xlabel('Step [k]')
ylabel('$H_{1,1}$ [s cm$^{-2}$]')

subplot(2,2,2)
plot(H(1:2:end,2),'k','linewidth',1.1)
grid('on')
xlim([0 n])
xlabel('Step [k]')
ylabel('$H_{1,2}$ [s cm$^{-2}$]')

subplot(2,2,3)
plot(H(2:2:end,1),'k','linewidth',1.1)
grid('on')
xlim([0 n])
xlabel('Step [k]')
ylabel('$H_{2,1}$ [s cm$^{-2}$]')

subplot(2,2,4)
plot(H(2:2:end,2),'k','linewidth',1.1)
ax = gca;
ax.YAxis.Exponent = -3;
grid('on')
xlim([0 n])
xlabel('Step [k]')
ylabel('$H_{2,2}$ [s cm$^{-2}$]')

sgtitle('Markov parameters of approximated system','FontSize',11)

set(subplot(2,2,1), 'Position', [.055, .55, .39, .35]);
set(subplot(2,2,2), 'Position', [.555, .55, .39, .35]);
set(subplot(2,2,3), 'Position', [.055, .085, .39, .35]);
set(subplot(2,2,4), 'Position', [.555, .085, .39, .35]);


