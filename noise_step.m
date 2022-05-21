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

%% Simulation Scenario - Noise response

% Simulation time frame
t_0 = 0;
t_f = 1300;
timeframe = t_0:t_f;
t_step = 0;

close all
for nl = [0.01 0.1 1] % Low, Medium, High noise
    
    pch = 1.1;
    [y_1_10,u_1_10,t_1_10] = noisesim(pch,1,u,d,p,xs,t_0,t_f,t_step,nl);
    [y_2_10,u_2_10,t_2_10] = noisesim(pch,2,u,d,p,xs,t_0,t_f,t_step,nl);
    
    pch = 1.25;
    [y_1_25,u_1_25,t_1_25] = noisesim(pch,1,u,d,p,xs,t_0,t_f,t_step,nl);
    [y_2_25,u_2_25,t_2_25] = noisesim(pch,2,u,d,p,xs,t_0,t_f,t_step,nl);
    
    pch = 1.5;
    [y_1_50,u_1_50,t_1_50] = noisesim(pch,1,u,d,p,xs,t_0,t_f,t_step,nl);
    [y_2_50,u_2_50,t_2_50] = noisesim(pch,2,u,d,p,xs,t_0,t_f,t_step,nl);
    
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
    
    sgtitle('Noisy response to percentage changes in input from steady state','FontSize',11)
    
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
    
    sgtitle('Normalized noisy response to percentage changes in input from steady state','FontSize',11)
    
    set(subplot(2,2,1), 'Position', [.055, .53, .39, .35]);
    set(subplot(2,2,2), 'Position', [.555, .53, .39, .35]);
    set(subplot(2,2,3), 'Position', [.055, .085, .39, .35]);
    set(subplot(2,2,4), 'Position', [.555, .085, .39, .35]);
    
end












