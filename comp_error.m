clear
clc

load('mpced.mat')
mpced = abs(ed);

load('pided.mat')
pided = abs(ed);

figure('Units','inches','Position',[0 0 6.693 4],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

subplot(2,1,1)
stairs(pided(:,1),'k','linewidth',1.3)
hold on
stairs(mpced(:,1),'Color', [0.5 0.5 0.5],'linewidth',1.3)
hold off
grid('on')
xlabel('Time [s]')
ylabel('$|e|$ [cm]')
xlim([0 250])
ylim([0 80])
legend('$\textrm{PID}$','$\textrm{MPC}$','interpreter','latex')

subplot(2,1,2)
stairs(pided(:,2),'k','linewidth',1.3)
hold on
stairs(mpced(:,2),'Color', [0.5 0.5 0.5],'linewidth',1.3)
hold off
grid('on')
xlabel('Time [s]')
ylabel('$|e|$ [cm]')
xlim([0 250])
ylim([0 80])
legend('$\textrm{PID}$','$\textrm{MPC}$','interpreter','latex')

sgtitle('Absolute errors between the reference and outputs','FontSize',11)

set(subplot(2,1,1), 'Position', [.065, .55, .90, .35]);
set(subplot(2,1,2), 'Position', [.065, .085, .90, .35]);