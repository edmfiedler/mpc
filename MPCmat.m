function [H,Mx0,Mr,Mu_1,Md,Aqp,Phi,Gamma_d] = MPCmat(A,B,C,Bv,Q,S,n)

Phi = zeros(n*size(C,1),size(C,2));
Phi_i = C;
for i = 1:size(C,1):n*size(C,1)
    Phi(i:i+size(C,1)-1,:) = Phi_i*A;
    Phi_i = Phi(i:i+size(C,1)-1,:);
end

% Obtain a vector of Markov Parameters
Gamma_vec = markovvec(A,B,C,n);
Gamma = zeros(size(Gamma_vec,1));
% Bulid Gamma Matrix
k = 0;
for i = 1:size(C,1):length(Gamma)
    Z = zeros(k*size(C,1),size(C,1));
    Gamma(:,i:i+size(C,1)-1) = [Z;Gamma_vec(1:end-size(Z,1),:)];
    k = k + 1;
end

% Obtain a vector of Markov Parameters
Gamma_vec = markovvec(A,Bv,C,n);
Gamma_d = zeros(size(Gamma_vec,1));
% Bulid Gamma Matrix
k = 0;
for i = 1:size(C,1):length(Gamma)
    Z = zeros(k*size(C,1),size(C,1));
    Gamma_d(:,i:i+size(C,1)-1) = [Z;Gamma_vec(1:end-size(Z,1),:)];
    k = k + 1;
end

Qz = zeros(length(Gamma));
k = 0;
for i = 1:length(Q):length(Gamma)
    Z_b = zeros(k*length(Q),length(Q));
    Z_a = zeros(length(Gamma)-size(Z_b,1)-length(Q),length(Q));
    Qz(:,i:i+length(Q)-1) = [Z_b; Q; Z_a];
    k = k + 1;
end

Hs = zeros(length(Gamma));
k = 0;
for i = 1:length(S):length(Gamma)
    if i == 1
        Hs_i = [2*S;-S];
        Hs(:,i:i+length(S)-1) = [Hs_i;zeros(length(Gamma)-size(Hs_i,1),length(S))];
    else
        if i == length(Gamma)-1
            Hs_i = [-S;S];
            Hs(:,i:i+length(S)-1) = [zeros(length(Gamma)-size(Hs_i,1),length(S));Hs_i];
        else
            Hs_i = [-S;2*S;-S];
            Z_b = zeros(k*length(S),length(S));
            Z_a = zeros(length(Gamma)-size(Z_b,1)-size(Hs_i,1),length(S));
            Hs(:,i:i+length(S)-1) = [Z_b;Hs_i;Z_a];
            k = k + 1;
        end
    end
end

Lam = zeros(length(Gamma));
k = 0;
for i = 1:length(Q):length(Gamma)
    if i == length(Gamma)-1
        Lam(:,i:i+length(Q)-1) = [zeros(length(Gamma)-length(Q),length(Q));eye(length(Q))];
    else
        Lam_i = [eye(length(Q)); -eye(length(Q))];
        Z_b = zeros(k*length(Q),length(Q));
        Z_a = zeros(length(Gamma)-size(Z_b,1)-size(Lam_i,1),length(Q));
        Lam(:,i:i+length(Q)-1) = [Z_b; Lam_i ; Z_a];
        k = k + 1;
    end
end

% Compute matrices
H = Gamma'*Qz*Gamma+Hs;
Mx0 = Gamma'*Qz*Phi;
Mr = -Gamma'*Qz;
Md = Gamma'*Qz*Gamma_d;
Mu_1 = -[S;zeros(length(Gamma)-length(S),length(S))];
Aqp = [Lam;Gamma];

end

