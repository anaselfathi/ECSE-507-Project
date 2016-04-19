function [ out ] = Betts1220ConstEq( X )
%BETTS6FUN Summary of this function goes here
%   Detailed explanation goes here

X = X(:);

if(length(X) ~= 16)
    out = NaN;
    return;
end

out(1) = 0.22*X(1) +0.20*X(2) +0.19*X(3) +0.25*X(4) +0.15*X(5) +0.11*X(6)+0.12*X(7) +0.13*X(8) + X(9) - 2.5;
out(2) = -1.46*X(1) - 1.3*X(3) + 1.82*X(4) - 1.15*X(5) + 0.8*X(7) + X(10) - 1.10;
out(3) = 1.29*X(1) - 0.89*X(2) - 1.16*X(5) - 0.96*X(6)- 0.49*X(8) + X(11) + 3.10;
out(4) = -1.10*X(1) - 1.06*X(2) + 0.95*X(3) - 0.54*X(4) - 1.78*X(6) - 0.41*X(7) + X(12) + 3.50;
out(5) = -1.43*X(4) + 1.51*X(5) + 0.59*X(6) - 0.33*X(7) - 0.43*X(8) + X(13) - 1.30;
out(6) = -1.72*X(2) - 0.33*X(3) + 1.62*X(5) + 1.24*X(6) + 0.21*X(7) - 0.26*X(8) + X(14) - 2.10;
out(7) = 1.12*X(1) + 0.31*X(4) + 1.12*X(7) - 0.36*X(9) + X(15) - 2.30;
out(8) = 0.45*X(2) + 0.26*X(3) - 1.10*X(4) + 0.58*X(5) - 1.03*X(7) + 0.1*X(8) + X(16) + 1.50;

out = -out(:);

end
