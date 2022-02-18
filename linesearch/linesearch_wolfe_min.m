function [x_min,x_it]=linesearch_wolfe_min(f,df,hf,x0,varargin)
% LINESEARCH_ARMIJO_MIN implements line-search minimization with a 
%   given funciton, gradient, Hessian, and initial guess.  Input arguments
%   are:
%
%   varargin:   optional arguments to set default parameters
%
%                  Parameter    Default     Description
%                   'r'          0.75        line parameter
%                   'eta'        0.50        wolfe parameter
%                   'mu'         0.50        optimization indicator
%                   'ep'         1E-5        desired accuracy
%                   'itmax'      20          maximum iterations
%                   'prt'        0           print details
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %initial parameters
    r     = 0.75;
    eta   = 0.50;
    mu    = 0.50;
    ep    = 1E-5;
    itmax = 50;
    prt   = 0;

    %parse variable input arguments
    if (~isempty(varargin))
        for c=1:length(varargin)
            switch varargin{c}
                case {'mu'}
                    mu=varargin{c+1};
                case {'r'}
                    r=varargin{c+1};
                case {'eta'}
                    eta=varargin{c+1};    
                case {'ep'}
                    ep=varargin{c+1};
                case {'itmax'}
                    itmax=varargin{c+1};
                case {'prt'}
                    prt=varargin{c+1};
            otherwise         
                if(~isnumeric(varargin{c}))
                    error(['Invalid optional argument, ', ...
                        varargin{c}]);
                end
            end
        end
    end
    
    %set initial guess value and make vector of x values per iteration
    x_size    = length(x0);
    x         = x0;
    x_it      = zeros(itmax+1,x_size);
    x_it(1,:) = x0;
    it        = 1;

    %check for errors before starting assignment and iteration
    if(mu<=0||mu>=1||r<=0||r>=1||eta<=0||eta>=1)
        error('Parameters must be defined such that 0 < mu, r, eta < 1');
    elseif(length(df(x))~=x_size||size(hf(x),1)~=x_size|| ...
            size(hf(x),2)~=x_size)
        error(['The initial guess, function, gradient, and Hessian' ...
            'must be of consistent dimension']);
    end
    
    %main loop
    while(norm(df(x))>ep&&(it<=itmax))

        %initialize line search parameter and find search direction
        alpha=1;
        p=-hf(x)\df(x);
        x_k=x+alpha*p;
        
        % Armijo condition
        while(abs(p'*df(x_k))>eta*abs(p'*df(x))&&alpha>ep)
            alpha = r*alpha;
            x_k=x+alpha*p;
        end
        
        % adjust x if the condition is met.
        if abs(p'*df(x_k))<=eta*abs(p'*df(x))
            x=x+alpha*p;
            x_it(it+1,:) = x';
            it=it+1;
        else 
            if prt; fprintf('alpha is too small\n'); end
            it=itmax+1;
        end

    end

    %final assignments
    x_min = x;
    x_it  = x_it(1:it,:);
    if it>itmax
        if prt; fprintf('reached the max iterations\n');end
    end
    end