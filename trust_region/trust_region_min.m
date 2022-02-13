function [x_min,x_it] = trust_region_min(f,df,hf,x0,varargin)
% TRUST_REGION_MIN uses trust-region-minimization on a function with
%   known gradient, Hessian, and initial guess.  Input arguments are
%   
%   f:          function of n variables in vector length n
%   df:         gradient in a 1-by-n size vector
%   hf:         Hessian in n-by-n vector
%
%   varargin:   optional arguments to set default parameters
%
%                  Parameter    Default     Description
%                   'trb'        1.0         trust region bound
%                   'mu'         0.50        optimization indicator
%                   'eta'        0.75        t-r upper criteria
%                   'ep'         1E-5        desired accuracy
%                   'itmax'      20          maximum iterations
%                   'prt'        0           print details
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %initial parameters
    mu      = 0.50;
    eta     = 0.75;
    trb     = 1.0;
    ep      = 1E-5;
    it      = 1;
    itmax   = 20;
    prt     = 0;
    
    %parse variable input arguments
    if (~isempty(varargin))
        for c=1:length(varargin)
            switch varargin{c}
                case {'mu'}
                    mu=varargin{c+1};
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
    x_size  = length(x0);
    x       = x0;
    x_it    = zeros(itmax+1,x_size);
    x_it(1,:) = x0;

    %check for errors before starting assignment and iteration
    if(mu<=0||mu>=1||eta<=mu||eta>=1)
        error('Parameters must be defined such that 0 < mu < eta < 1');
    elseif(length(df(x))~=x_size||size(hf(x),1)~=x_size|| ...
            size(hf(x),2)~=x_size)
        error(['The initial guess, function, gradient, and Hessian' ...
            'must be of consistent dimension']);
    end
    
    %functions used to assess predicted change
    phi=@(x,p) f(x)+df(x)'*p+0.5*p'*hf(x)*p;
    rho =@(x,p) (f(x)-f(x+p))/(f(x)-phi(x,p));
    
    %main loop
    while(norm(df(x))>ep && it<itmax)
    
        %determine predicted change p and rescale if greater than threshold
        p=-hf(x)\df(x);
        if trb<norm(p)
            g=@(l) norm((hf(x)+l*eye(x_size))\df(x))-trb;
            lambda = fzero(g,0.9);
            p=-(hf(x)+lambda*eye(x_size))\df(x);
        end
        
        %update x if p is significant
        pk = rho(x,p);
        if pk > mu
            x=x+p;
        end
        
        %update threshold if needed
        if pk <=mu
            trb = 0.5*trb;
        elseif pk>=eta
            trb = 2*trb;
        end

        %
        x_it(it+1,:) = x';
        it       = it+1;
    end 
    %assign final minimum return value
    x_min = x;
    x_it  = x_it(1:it,:);
    
    %print out options if flag set
    if(prt)
        %if the maximum number of iterations is met before convergence
        if(it>=itmax)
            fprintf('no convergence within %3.0f iterations\n',itmax)
            fprintf('final divergence: %4.8f\n',norm(df(x)))
    
        %if the function converged within tolerance    
        else
            fprintf('iterations: %3.0f\n',it-1)
        end
    end
end