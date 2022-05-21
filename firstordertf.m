function [tf_1,tf_2,K,tau] = firstordertf(y,u,t,hs,t_step)
check = u(1,:) == u(end,:);
st = find(0 == check);

K = [(y(end,1)-hs(1))/(u(end,st)-u(1,st)) (y(end,2)-hs(2))/(u(end,st)-u(1,st))];

ind = find(y(:,1)<(y(end,1)-hs(1))*0.63+hs(1));
tau_1 = t(ind(end))-t_step;

ind = find(y(:,2)<(y(end,2)-hs(2))*0.63+hs(2));
tau_2 = t(ind(end))-t_step;

tau = [tau_1 tau_2];

tf_1 = tf(K(1),[tau_1 1]);
tf_2 = tf(K(2),[tau_2 1]);
end

