function [ out ] = Betts1110Const( X )
%Betts8Const Summary of this function goes here
%   Detailed explanation goes here
X = X(:);
if(length(X) ~= 2)
    out = NaN;
    return;
end

out(1) = X(1)*X(1) + X(2)*X(2) - 25;
out(2) = X(1)*X(2) - 9;

out = out(:);

end

