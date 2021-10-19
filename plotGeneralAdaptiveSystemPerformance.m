function plotAdaptiveSystemPerformance(t_ode,xr,xd,x,V1,V2,V3,thx,thr,thx_star,thr_star,tplotrange,plotEveryNPoints)
SubpFS = 8.5;
SFigFS = 11;
DFigFS = 12;

pp = plotEveryNPoints;

figure('Renderer', 'painters', 'Position', [0 10 600 800]); clf
subplot(5,1,1)
    set(gca,'FontSize',SubpFS); hold on
    plot(t_ode(1:pp:end),xr(1,1:pp:end),'r',t_ode(1:pp:end),xd(1,1:pp:end),'k--',t_ode(1:pp:end),x(1,1:pp:end),'b')
    legend('x_r=\omega_r','\omega_d','x=\omega','Location','northeast')
    xlabel('Time (s)')
    ylabel('Angular Velocity')
    xlim(tplotrange)
    hold off
    box on 
subplot(5,1,2)
    set(gca,'FontSize',SubpFS); hold on
    plot(t_ode(1:pp:end),x(:,1:pp:end)-xr(:,1:pp:end),'r')
    ylabel('Error')
    xlabel('Time (s)')
    legend('e=x-x_r','Location','southeast')
    ylim([min(min(x-xr))-.1 max(max(x-xr))+.1])
    xlim(tplotrange)
    grid on
    hold off
    box on
subplot(5,1,3)
    set(gca,'FontSize',SubpFS); hold on
    plot(t_ode(1:pp:end),thx(:,1:pp:end),'r-',[0,t_ode(end)],thx_star*[1 1],'r--')
    legend('\theta_x','','','','\theta_x*','','','','Location','southeast')
    xlabel('Time (s)')
    ylabel('Adaptive Parameters')
    %ylim([min(min(thx),min(thr))-.1 max(max(thx),max(thr))+.03])
    xlim(tplotrange)
    hold off
    box on 
subplot(5,1,4)
    set(gca,'FontSize',SubpFS); hold on
    plot(t_ode(1:pp:end),thr(:,1:pp:end),'b-')
    for i = 1:length(thr_star)
        plot([0,t_ode(end)],thr_star(i)*[1 1],'b--')
    end
    legend('\theta_B','','\theta_B*','','Location','southeast')
    xlabel('Time (s)')
    ylabel('Adaptive Parameters')
    %ylim([min(min(thx),min(thr))-.1 max(max(thx),max(thr))+.03])
    xlim(tplotrange)
    hold off
    box on 
subplot(5,1,5)
    set(gca,'FontSize',SubpFS); hold on
    plot(t_ode(1:pp:end),V1(1:pp:end),'g',t_ode(1:pp:end),V2(1:pp:end),'r--',t_ode(1:pp:end),V3(1:pp:end),'b:')
    legend('V_1','V_2','V_3')
    ylabel('Lyapunov Fcn. V(e,\theta_x,\theta_r)')
    xlabel('Time (s)')
    xlim(tplotrange)
    hold off
    box on 

end