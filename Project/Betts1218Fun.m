function [ out ] = Betts1218Fun( X )
%BETTS6FUN Summary of this function goes here
%   Detailed explanation goes here

X = X(:);
if(length(X) ~= 5)
    out = NaN;
    return;
end

out = 5.3578547 * X(3)^2 + 0.8356891*X(1)*X(5) + 37.293239*X(1) - 40792.141;

end

