function [y,ud,Td] = detsim(pch,u_i,u,d,p,xs,t_0,t_f,t_step)
% Simulation time frame
timeframe = t_0:t_f;

y = []; Td = []; ud = []; x0 = xs;
for time = timeframe
    [T,X] = ode15s(@ModifiedFourTankSystem,[time time+1],x0,[],u,d,p);
    x0 = X(end,:);
    Td = [Td;T];
    y = [y;(1/(p(5)*p(12))).*X(:,1) (1/(p(6)*p(12))).*X(:,2)];
    z = y;
    ud = [ud; u' d'];
    
    if time == t_step
        if u_i == 1
            u = [u(1)*pch; u(2)];
        end
        
        if u_i == 2
            u = [u(1); u(2)*pch];
        end
    end
end

end

