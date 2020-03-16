function [X,Wm,Wc] = sigmaPoints(x_kk, p_kk, alpha, beta, ki)
%SIGMAPOINTS Choose sigma points for UKF
L = numel(x_kk);
lambda=alpha^2*(L+ki)-L;
Wm=[lambda/(lambda+L) 0.5/(lambda+L)*ones(1,2*L)];
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);
A = chol(p_kk);
X = ones(3,2*L+1);
X(:,1) = x_kk;

for i=2:L+1
    X(:,i) = x_kk + sqrt(L+lambda)*A(i-1,:)';
    X(:,i+L) = x_kk - sqrt(L+lambda)*A(i-1,:)';
end
end

