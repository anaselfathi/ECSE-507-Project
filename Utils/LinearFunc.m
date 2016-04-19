function out = LinearFunc( X, a, b)
%LinearFunc Summary of this function goes here
%   a^T X = b

X = X(:);
a = a(:);
b = b(:);

if(sum(size(X) ~= size(a)) > 0)
    error('myfuns:LinearFunc:WrongInput', ...
        'X and a should have the same size');
end

out = a'*X + b;

end

