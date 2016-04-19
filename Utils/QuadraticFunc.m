function out = QuadraticFunc( X, A, b, c)
%QUADRATIC Summary of this function goes here
%   X^T A X + b^T X + c

X = X(:);
b = b(:);
c = c(:);

if(size(A,1) ~= size(A,2))
    error('myfuns:QuadraticFunc:WrongInput',...
         'Q is squared');
end

if(size(A,1) ~= size(X,1))
    error('myfuns:QuadraticFunc:WrongInput',...
        'X and Q should have the same size');
end

if(length(b) ~= length(X))
    error('myfuns:QuadraticFunc:WrongInput',...
        'X and b should have the same size');
end

out = X' * A * X + b' * X + c;

end

