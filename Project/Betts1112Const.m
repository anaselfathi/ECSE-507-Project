function [ out ] = Betts1112Const( X )
%Betts8Const Summary of this function goes here
%   Detailed explanation goes here
X = X(:);
if(length(X) ~= 2)
    out = NaN;
    return;
end

out = 10*(X(2) - X(1))^2;

end

