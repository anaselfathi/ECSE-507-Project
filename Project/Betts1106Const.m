function [ out ] = Betts1106Const( X )
%BETTS6FUN Summary of this function goes here
%   Detailed explanation goes here

X = X(:);

if(length(X) ~= 7)
    out(1) = NaN;
    out(2) = NaN;
    return;
end

if(sum(X > pi/2) ~= 0)
    out(1) = inf;
    out(2) = inf;
    return;
end

a = [0;50;50;75;75;75;100;100];
t = [0;25;50;100;150;200;290;380];
g = 32;

p = zeros(1,8);
pd = zeros(1,8);
q = zeros(1,8);
qd = zeros(1,8);

for i = 2:1:8
    p(i) = 0.5 * a(i) * (t(i) - t(i-1))^2 * cos(X(i-1)) + pd(i-1) * (t(i) - t(i-1)) + p(i-1);
    q(i) = 0.5 * (t(i) - t(i-1))^2 * (a(i) * sin(X(i-1)) - g) + qd(i-1) * (t(i) - t(i-1)) + q(i-1);
    pd(i) = a(i)*(t(i) - t(i-1)) * cos(X(i-1)) + pd(i-1);
    qd(i) = (t(i) - t(i-1)) * (a(i) * sin(X(i-1)) - g) + qd(i-1);
end

out(1) = q(8) - 100000;
out(2) = qd(8) - 1000;

out = out(:);

end

