function [ Xopt, X, V, n ] = SteepestDescent(f, df, X0, N)
%STEEPESTDESCENT_OPT Summary of this function goes here
%   Detailed explanation goes here

X0 = X0(:);

X = NaN*ones(length(X0), N);
V = NaN*ones(1, N);

gamma = 1.5;
mu = 0.8;

n = 1;
X(:,n) = X0;
V(n) = f(X(:,n));
dvX = df(X(:,n));

while n < N
    if(norm(dvX) < 1e-8)
        break;
    end
    s = -dvX;
    w = armijoLineSearch(f, dvX, X(:,n), s, gamma, mu, 100);
   
    if(norm(w) < 1e-8)
        break;
    end
        
    X(:,n+1) = X(:,n) + w*s;
    n = n + 1;
    
    V(n) = f(X(:,n));
    dvX = df(X(:,n));
end

Xopt = X(:,n-1);
if(sum(sum(isnan(X))))
    X(isnan(X)) = [];
    X = reshape(X, length(X0), length(X)/length(X0));
    V(isnan(V)) = [];
end

end

