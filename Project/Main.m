clc,
close all,
addpath('../Utils');
format compact
format long
%% simple example
f = @(x)(LinearFunc(x, [1;1], 0));
df = @(x)([1;1]);
h = @(x)(QuadraticFunc(x, eye(2), [0;0], -1));
dh =  @(x)(2*x);
g = @(x)(-x);
dg = @(x)([-1 0;0 -1]);

%% f & h
[ExBFGS.X, ExBFGS.steps] = AugmentedLagrangian(f,h,[0;0],0,...
    'epsState', 1e-5,...
    'epsCost', 1e-10,...
    'Method', 'bfgs');

%% f & g
[ExDFP.X, ExDFP.steps] = AugmentedLagrangian(f,[],g,[0;0],[],[0;0],...
    'epsState', 1e-5, ...
    'epsCost', 1e-6,...
    'Method', 'dfp');

%% f & h & g
[ExSymRank.X, ExSymRank.steps] = AugmentedLagrangian(f,h,g,[0;0],0,[0;0],...
    'epsState', 1e-5,...
    'epsCost', 1e-6,...
    'Method', 'symRank');

%% f & h & g & df & dh & dg
[ExSymRankDiff.X, ExSymRankDiff.steps] = AugmentedLagrangian(f,h,g,[0;0],0,[0;0],...
    'epsState', 1e-4, ...
    'epsCost', 1e-8,...
    'Method', 'symRank',...
    'fdiff',df, ...
    'gdiff',dg,...
    'hdiff',dh);

%% Betts 11.6
[Betts1106.X, Betts1106.Steps] = AugmentedLagrangian(@Betts1106Fun, @Betts1106Const, 0.5*ones(7,1), [0;0],...
    'Method', 'dfp',...
    'penalty', 75,...
    'MaxIter', 25,...
    'epsState', 1e-9,...
    'epsCost', 1e-5);

save('BettsOptim','Betts1106');
%% Betts 11.8
[Betts1108.X, Betts1108.steps] = AugmentedLagrangian(@Betts1108Fun, @Betts1108Const, [0.7;0.2;0.1], 0,...
    'MaxIter', 25,...
    'epsState', 1e-5,...
    'epsCost', 1e-5,...
    'Method', 'bfgs',...
    'penalty', 1e4);

save('BettsOptim','Betts1108','-append');
%% Betts 11.10
[Betts1110.X, Betts1110.steps] = AugmentedLagrangian(@Betts1110Fun, @Betts1110Const, [2;2], [0;0],...
    'MaxIter', 25,...
    'epsState', 1e-10,...
    'epsCost', 1e-20,...
    'Method', 'bfgs',...
    'penalty', 1e15, 'verbose', 3);

save('BettsOptim','Betts1110','-append');
%% Betts 11.12
[Betts1112.X, Betts1112.steps] = AugmentedLagrangian(@Betts1112Fun, @Betts1112Const, [-1.2;1.0], 0,...
    'penalty', 5,...
    'MaxIter', 25,...
    'epsState', 1e-10,...
    'epsCost', 1e-20,...
    'Method', 'bfgs');

save('BettsOptim','Betts1112','-append');
%% Betts 12.1
[Betts1201.X, Betts1201.steps] = AugmentedLagrangian(@Betts1201Fun,[],@Betts1201Const,[2;2;2;2;2],[],zeros(10,1),...
    'MaxIter', 25,...
    'penalty', 50,...
    'epsState', 1e-8,...
    'epsCost', 1e-8,...
    'Method', 'dfp');

save('BettsOptim','Betts1201','-append');
%% Betts 12.8
[Betts1208.X, Betts1208.steps] = AugmentedLagrangian(@Betts1208Fun,[],@Betts1208Const,[-2;1],[],zeros(5,1),...
    'MaxIter', 25,...
    'penalty', 1,...
    'epsState', 1e-6,...
    'epsCost', 1e-10,...
    'Method', 'dfp');

save('BettsOptim','Betts1208','-append');
%% Betts 12.18
[Betts1218.X, Betts1218.steps] = AugmentedLagrangian(@Betts1218Fun,[],@Betts1218Const,[78.62;33.44;31.07;44.18;35.32],[],zeros(16,1),...
    'MaxIter', 25,...
    'penalty', 100,...
    'epsState', 1e-6,...
    'epsCost', 1e-8,...
    'Method', 'bfgs');

save('BettsOptim','Betts1218','-append');
%% Betts 12.20
[Betts1220.X, Betts1220.steps] = AugmentedLagrangian(@Betts1220Fun, @Betts1220ConstEq, @Betts1220ConstIneq, zeros(16,1), zeros(8,1), zeros(32,1),...
    'MaxIter', 15,...
    'penalty', 100,...
    'epsState', 1e-6,...
    'epsCost', 1e-4,...
    'Method', 'dfp');

save('BettsOptim','Betts1220','-append');
%% Betts 12.22
[Betts1222.X, Betts1222.steps] = AugmentedLagrangian(@Betts1222Fun, [], @Betts1222Const, [10;10;10], [], zeros(8,1),...
    'MaxIter', 25,...
    'penalty', 1000,...
    'epsState', 1e-5,...
    'epsCost', 1e-8,...
    'Method', 'bfgs');

save('BettsOptim','Betts1222','-append');
%% Betts 12.27
[Betts1227.X, Betts1227.steps] = AugmentedLagrangian(@Betts1227Fun, [], @Betts1227Const, [2.52;2;37.5;9.25;6.8], [], zeros(15,1), ...
    'fdiff', @Betts1227FunDiff, ...
    'gdiff', @Betts1227ConstDiff,...
    'MaxIter', 100,...
    'penalty', 1e5,...
    'penaltyMulti', 2.5,...
    'epsState', 1e-15,...
    'epsCost', 1e-22,...
    'Method', 'bfgs');

save('BettsOptim','Betts1227','-append');
%%

% options = optimoptions(@fmincon,'Algorithm','sqp','display','testing');
% [x,fval] = fmincon(@Betts1227Fun,[2.52;2;37.5;9.25;6.8],[],[],[],[],[],[],...
%    @(x) deal(Betts1227Const(x), []),options);