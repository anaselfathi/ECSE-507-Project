function diff = numDiff(f, X, eps)
%NUMDIFF Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 3)
    eps = 1e-4;
end

h = abs(X) * eps;
h(h < eps*1e-1) = eps*1e-1;

X = X(:);
h = h(:);

if(length(X) ~= length(h))
    error('myfuns:numDiff:WrongInput', ...
        'step h should have the same size as X');
end

diff = NaN * ones(1, length(X));

I = eye(length(X));

for i = 1:length(X)
    diff(i) = (f(X + h(i)*I(:,i)) - f(X))/h(i);
end

diff = diff(:);

end

