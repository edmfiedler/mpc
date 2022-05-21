function [xh,Pk] = statestimator(xh,Pk,A,B,Bv,C,u,y,R,Q)
    Rek = C*Pk*C'+R;
    Kf = Pk*C'*inv(Rek);
    ek = y - C*xh;
    xh = xh + Kf*ek;
    Pk = Pk - Kf*Rek*Kf';
    
    xh = A*xh+B*u;
    Pk = A*Pk*A'+Bv*Q*Bv';
end

