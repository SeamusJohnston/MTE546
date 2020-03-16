function [X,W] = sigmaPoints(x_kk, p_kk, W0)
%SIGMAPOINTS Choose sigma points for UKF
n = numel(x_kk);
A = chol(p_kk);
X = ones(3,2*n+1);
X(:,1) = x_kk;

for i=2:n+1
    X(:,i) = x_kk + sqrt(n/(1-W0))*A(:,i-1);
    X(:,i+n) = x_kk - sqrt(n/(1-W0))*A(:,i-1);
end

W = (1-W0)/(2*n) * ones(1,2*n);
W = [W0 W];
end

