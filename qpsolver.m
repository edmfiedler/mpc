function [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit)
    A = [A;-A];
    b = [bu;-bl];   
    [x,info] = quadprog(H,g,A,b,[],[],l,u,xinit);
end

