function [ Xopt, Steps ] = Secant(f, df, X0, varargin)
%Secant: Quasi-Newton algorithm for unconstrained Optimization
%   min f
%
% f: objective function
%
% Usage:
% - Secant(f,[],X0, ...)
% - Secant(f,df,X0, ...)
%
% - Optional parameter:
%   + maxiter: Max iterations to run
%   + method: Constrained optimization method,
%          supported: bfgs, dfp, symmetric rank
%   + numdiff: Force numeric differentiation
%   + epsstate: Error on X
%   + epsCost: Error on Constraint function
%   + verbose: Verbose level: 0, 1, 2.
%
% author: anas.elfathi@mail.mcgill.ca - 2016
tic;

% default option
N = 100;
method = 0;
if(isempty(df))
    numdiff = 1;
else
    numdiff = 0;
end
eps1 = 1e-5;
eps2 = 1e-10;
verbose = 0;
methodStr = {'BFGS', 'DFP', 'Symmetric Rank'};

% options
for nVar = 1:2:length(varargin)
    switch(lower(varargin{nVar}))
        case 'maxiter'
            N = varargin{nVar+1};
        case 'method'
            switch(lower(varargin{nVar+1}))
                case 'bfgs'
                    method = 0;
                case 'dfp'
                    method = 1;
                case 'symrank'
                    method = 2;
                otherwise
                    error('myfuns:Secant:WrongInput', ...
                        'Unkown method');
            end
        case 'numdiff'
            switch(lower(varargin{nVar+1}))
                case {'yes', 'y', 1}
                    numdiff = 1;
                otherwise
                    numdiff = 0;
            end
        case 'epsstate'
            eps1 = varargin{nVar+1};
        case 'epscost'
            eps2 = varargin{nVar+1};
        case 'verbose'
            switch(lower(varargin{nVar+1}))
                case {'testing', 'test', 't', 2}
                    verbose = 2;
                case {'notify', 'yes', 'y', 1}
                    verbose = 1;
                otherwise
                    verbose = 0;
            end
        otherwise
            error('myfuns:Secant:WrongInput', ...
                'Unkown option %s', varargin{nVar});
    end
end

X0 = X0(:);

Steps.X = NaN*ones(length(X0), N);
Steps.f = NaN*ones(1, N);
Steps.df = NaN*ones(length(X0), N);

tet = 75*pi/180;
Gamma = 1.5:0.05:2.05;
Mu = 0.5:0.05:0.95;

n = 1;
Steps.X(:,n) = X0;
Steps.f(n) = f(X0);
if(~numdiff)
    Steps.df(:,n) = df(X0);
else
    Steps.df(:,n) = numDiff(f, X0, eps1);
end
Steps.H{n} = ((eps1)^2/eps2)*eye(length(X0));

LineVerbose{n} = sprintf(               '+------+-----------------------------Secant-%15s, Max iterations %7d: -----------------------------------+\\n', methodStr{method+1}, N);
LineVerbose{n} = strcat(LineVerbose{n}, '| Iter |                  State                   |   Function    |                 Gradient                 |    Step   |\n');
LineVerbose{n} = strcat(LineVerbose{n}, '+------+------------------------------------------+---------------+------------------------------------------+-----------+\n');
LineVerbose{n} = strcat(LineVerbose{n}, sprintf('|%06d| %40s | %+8.6e | %40s | --------- | \\n', n, Vector2String(Steps.X(:,n)), Steps.f(n), Vector2String(Steps.df(:,n))));

if(verbose > 0)
    fprintf(LineVerbose{n});
end

while n < N
    s = -Steps.H{n}*Steps.df(:,n);
    
    if(norm(s)*norm(Steps.df(:,n)) * cos(tet) > s' * (-Steps.df(:,n)))
        Steps.H{n} = ((eps1)^2/eps2)*eye(length(X0));
        s = -Steps.H{n}*Steps.df(:,n);
    end
    
    % line search best w, with resolution 1e-2
    k = 1;
    w = NaN * ones(1, length(Gamma) * length(Mu));
    vw = NaN * ones(1, length(Gamma) * length(Mu));
    for g = Gamma
        for m = Mu
            w(k) = ArmijoLineSearch(f, Steps.df(:,n), Steps.X(:,n), s, g, m, min(N*10, 1e2));
            vw(k) = f(Steps.X(:,n) + w(k)*s);
            k = k+1;
        end
    end
    [~, wIdx] = min(vw);
    wMin = w(wIdx);
    
    dx = wMin * s;
         
    if(norm(abs(dx)) < eps1)
        break;
    end
    
    Steps.X(:,n+1) = Steps.X(:,n) + dx;
    if(~numdiff)
        Steps.df(:,n+1) = df(Steps.X(:,n+1));
    else
        Steps.df(:,n+1) = numDiff(f, Steps.X(:,n+1), eps1);
    end
    
    if(max(abs(Steps.df(:,n+1))) < eps2)
        break;
    end
    
    if(mod(n, length(X0)) == 0)
        Steps.H{n+1} = ((eps1)^2/eps2)*eye(length(X0));
    else
        dg = Steps.df(:,n+1) - Steps.df(:,n);
        switch(method)
            case 0 %BFGS
                Steps.H{n+1} = Steps.H{n} + (1 + dg'*(Steps.H{n} * dg)/(dx' * dg)) * dx * dx'/(dx' * dg) - (dx*(Steps.H{n} * dg)' + (Steps.H{n} * dg)*dx')/(dx' * dg);
            case 1 %DFP
                Steps.H{n+1} = Steps.H{n} + (dx * dx') / (dx' * dg) - ((Steps.H{n} * dg) * (Steps.H{n} * dg)') / (dg' * Steps.H{n} * dg);
            case 2 % Symmetric Rank
                if(max(abs(dg)) < eps2)
                    Steps.H{n+1} = Steps.H{n};
                else
                    Steps.H{n+1} = Steps.H{n} + ((dx - Steps.H{n}*dg)*(dx - Steps.H{n}*dg)')/((dx - Steps.H{n}*dg)'*(dg));
                end
                
        end
    end
       
    n = n + 1;
    
    Steps.f(n) = f(Steps.X(:,n));
    
    LineVerbose{n} = sprintf('|%06d| %40s | %+8.6e | %40s | %+5.2e |\n', n, Vector2String(Steps.X(:,n)), Steps.f(n), Vector2String(Steps.df(:,n)), wMin);
    if(verbose > 1)
        fprintf(LineVerbose{n});
    end
    
    if(isnan(Steps.f(n)))
        error('myfuns:Secant:NaN', ...
            'Error in computations');
    end
    
    if(~isfinite(Steps.f(n)))
        error('myfuns:Secant:Infty', ...
            'Error in computations');
    end
end

if(verbose > 0)
    if(verbose == 1 && n > 1)
        fprintf(LineVerbose{n});
    end
    LineVerbose{n+1} = '+------+------------------------------------------+---------------+------------------------------------------+-----------+\n';
    fprintf(LineVerbose{n+1});
end

Xopt = Steps.X(:,n);
if(sum(sum(isnan(Steps.X))))
    Steps.X(isnan(Steps.X)) = [];
    Steps.X = reshape(Steps.X, length(X0), length(Steps.X)/length(X0));
    Steps.f(isnan(Steps.f)) = [];
    Steps.df(isnan(Steps.df)) = [];
    Steps.df = reshape(Steps.df, length(X0), length(Steps.df)/length(X0));
end
Steps.NumOfIterations = n;
Steps.Time = toc;

Steps.Verbose = sprintf(cell2mat(LineVerbose));

end

