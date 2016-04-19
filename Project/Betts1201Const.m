function [ out ] = Betts1201Const( X )
%Betts8Const Summary of this function goes here
%   Detailed explanation goes here
X = X(:);
if(length(X) ~= 5)
    out = NaN;
    return;
end

out = [-X;X-[1;2;3;4;5]];

end

