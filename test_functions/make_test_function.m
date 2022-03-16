function [f,g,h] = make_test_function(varargin)
% MAKE_TEST_FUNCTION returns handles to a specificed test function, as well
%   as its gradient and Hessian.
%
%   varargin:   specifies the test function desired 
%
%                  Function Name              Global Minimum 
%                   'ackley'                    f(0,0)       = 0
%                   'beale'                     f(3,0.5)     = 0
%                   'booth'                     f(1,3)       = 0
%                   'bulkinN6'                  f(-10,1)     = 0
%                   'goldstein-price'           f(0,-1)      = 3
%                   'himmelblau'                f(3,2)       = 0
%                                              ~f(-2.8,3.1)  = 0
%                                              ~f(-3.8,-3.3) = 0
%                                              ~f(3.6,-1.8)  = 0
%                   'holderTable'              ~f(8.1,9.7)   = -19.20
%                   'rastrigin2'                f(0,0)       = 0
%                   'rastrigin4'                f(0,0,0,0)   = 0
%                   'rosenbrock'                f(1,1)       = 0
%                   'sphere2'                   f(0,0)       = 0
%                   'sphere4'                   f(0,0,0,0)   = 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %parse variable input arguments
    if (~isempty(varargin))
        if(length(varargin)>2)
            error(['Use one argument for function name',varargin]);
        elseif(length(varargin)==2)
            if(varargin{2}=='v')
                vecflag = true;
            else
                error(['Second argument v specifies vectorization',...
                    varargin]);
            end
        else
            vecflag = false;
        end


        switch varargin{1}
            case {'ackley'}
                if vecflag
                    f=@ackley_funv;
                    g=@ackley_gfunv;
                    h=@ackley_hfunv;

                else
                    f=@ackley_fun;
                    g=@ackley_gfun;
                    h=@ackley_hfun;
                end


            case {'beale'}
                if vecflag
                    f=@beale_funv;
                    g=@beale_gfunv;
                    h=@beale_hfunv;

                else
                    f=@beale_fun;
                    g=@beale_gfun;
                    h=@beale_hfun; 
                end

            case {'booth'}
                if vecflag
                    f=@booth_funv;
                    g=@booth_gfunv;
                    h=@booth_hfunv;

                else
                    f=@booth_fun;
                    g=@booth_gfun;
                    h=@booth_hfun; 
                end

            case {'bulkinN6'}
                if vecflag
                    f=@bulkinN6_funv;
                    g=@bulkinN6_gfunv;
                    h=@bulkinN6_hfunv;

                else                
                    f=@bulkinN6_fun;
                    g=@bulkinN6_gfun;
                    h=@bulkinN6_hfun;    
                end

            case {'goldstein-price'}
                if vecflag
                    f=@gpf_funv;
                    g=@gpf_gfunv;
                    h=@gpf_hfunv;

                else                
                    f=@gpf_fun;
                    g=@gpf_gfun;
                    h=@gpf_hfun;
                end

            case {'himmelblau'}
                if vecflag
                    f=@himmelblau_funv;
                    g=@himmelblau_gfunv;
                    h=@himmelblau_hfunv;

                else
                    f=@himmelblau_fun;
                    g=@himmelblau_gfun;
                    h=@himmelblau_hfun;
                end

            case {'holdertable'}
                if vecflag
                    f=@holdertable_funv;
                    g=@holdertable_gfunv;
                    h=@holdertable_hfunv;

                else
                    f=@holdertable_fun;
                    g=@holdertable_gfun;
                    h=@holdertable_hfun; 
                end

            case {'rastrigin2'}
                if vecflag
                    f=@rastrigin2_funv;
                    g=@rastrigin2_gfunv;
                    h=@rastrigin2_hfunv;

                else
                    f=@rastrigin2_fun;
                    g=@rastrigin2_gfun;
                    h=@rastrigin2_hfun;
                end
               
            case {'rastrigin4'}
                if vecflag
                    f=@rastrigin4_funv;
                    g=@rastrigin4_gfunv;
                    h=@rastrigin4_hfunv;

                else
                    f=@rastrigin4_fun;
                    g=@rastrigin4_gfun;
                    h=@rastrigin4_hfun;
                end

            case {'rosenbrock'}
                if vecflag
                    f=@rosenbrock_funv;
                    g=@rosenbrock_gfunv;
                    h=@rosenbrock_hfunv;
                else
                    f=@rosenbrock_fun;
                    g=@rosenbrock_gfun;
                    h=@rosenbrock_hfun;
                end
                   
            case {'sphere2'}
                if vecflag
                    f=@sphere2_funv;
                    g=@sphere2_gfunv;
                    h=@sphere2_hfunv;
                else
                    f=@sphere2_fun;
                    g=@sphere2_gfun;
                    h=@sphere2_hfun;
                end

            case {'sphere4'}
                if vecflag
                    f=@sphere4_funv;
                    g=@sphere4_gfunv;
                    h=@sphere4_hfunv;
                else
                    f=@sphere4_fun;
                    g=@sphere4_gfun;
                    h=@sphere4_hfun;  
                end    

        otherwise         
            if(~isnumeric(varargin))
                error(['Invalid optional argument, ', ...
                    varargin]);
            end
        end
    end
end


% BEALE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z=beale_fun(x,y)
    z=(x*y^2-x+9/4)^2+(x*y^3-x+21/8)^2+(x*y-x+3/2)^2;
end

function z=beale_funv(x,y)
    z=(x.*y.^2-x+9./4).^2+(x.*y.^3-x+21./8).^2+(x.*y-x+3./2).^2;
end

function z=beale_gfun(x,y)
    z = [2*(y^2-1)*(x*y^2-x+9/4)+2*(y^3-1)*(x*y^3-x+21/8)+2*(y-1)*...
         (x*y-x+3/2);
          2*x*(x*y-x+3/2)+4*x*y*(x*y^2-x+9/4)+6*x*y^2*(x*y^3-x+21/8)];
end

function z=beale_gfunv(x,y)
    z = [2.*(y.^2-1).*(x.*y.^2-x+9./4)+2.*(y.^3-1).*(x.*y.^3-x+21./8)+...
        2.*(y-1).*(x.*y-x+3./2);
          2.*x.*(x.*y-x+3./2)+4.*x.*y.*(x.*y.^2-x+9./4)+6.*x.*y.^2.*...
          (x.*y.^3-x+21./8)];
end

function z=beale_hfun(x,y)
    z = [2*(y-1)^2+2*(y^2-1)^2+2*(y^3-1)^2,...
        6*y^2*(x*y^3-x+21/8)-2*x+2*x*y+2*x*(y-1)+4*y*(x*y^2-x+9/4)+4*x*...
        y*(y^2-1)+6*x*y^2*(y^3-1)+3;
         6*y^2*(x*y^3-x+21/8)-2*x+2*x*y+2*x*(y-1)+4*y*(x*y^2-x+9/4)+4*...
         x*y*(y^2-1)+6*x*y^2*(y^3-1)+3,...
         8*x^2*y^2+18*x^2*y^4+4*x*(x*y^2-x+9/4)+2*x^2+12*x*y*...
         (x*y^3-x+21/8)];
end

function z=beale_hfunv(x,y)
    z = [2.*(y-1).^2+2.*(y.^2-1).^2+2.*(y.^3-1).^2,...
        6.*y.^2.*(x.*y.^3-x+21./8)-2.*x+2.*x.*y+2.*x.*(y-1)+4.*y.*(x.*...
        y.^2-x+9./4)+4.*x.*y.*(y.^2-1)+6.*x.*y.^2.*(y.^3-1)+3;
         6.*y.^2.*(x.*y.^3-x+21./8)-2.*x+2.*x.*y+2.*x.*(y-1)+4.*y.*(x.*...
         y.^2-x+9./4)+4.*x.*y.*(y.^2-1)+6.*x.*y.^2.*(y.^3-1)+3,...
         8.*x.^2.*y.^2+18.*x.^2.*y.^4+4.*x.*(x.*y.^2-x+9./4)+2.*x.^2+...
         12.*x.*y.*(x.*y.^3-x+21./8)];
    z = reshape_hessian(z);
end


% GOLDSTEIN_PRICE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z=gpf_fun(x,y)
    z = ((x+y+1)^2*(3*x^2+6*x*y-14*x+3*y^2-14*y+19)+1)*((2*x-3*y)^2*...
        (12*x^2-36*x*y-32*x+27*y^2+48*y+18)+30);
 
end

function z=gpf_funv(x,y)
    z = ((x+y+1).^2.*(3.*x.^2+6.*x.*y-14.*x+3.*y.^2-14.*y+19)+1).*((2.*...
        x-3.*y).^2.*(12.*x.^2-36.*x.*y-32.*x+27.*y.^2+48.*y+18)+30);
 
end 

function z=gpf_gfun(x,y)
    z = [((6*x+6*y-14)*(x+y+1)^2+(2*x+2*y+2)*(3*x^2+6*x*y-14*x+3*y^2-14*...
         y+19))*((2*x-3*y)^2*(12*x^2-36*x*y-32*x+27*y^2+48*y+18)+30)+...
         ((x+y+1)^2*(3*x^2+6*x*y-14*x+3*y^2-14*y+19)+1)*((8*x-12*y)*(12*...
         x^2-36*x*y-32*x+27*y^2+48*y+18)-(2*x-3*y)^2*(36*y-24*x+32));
         ((6*x+6*y-14)*(x+y+1)^2+(2*x+2*y+2)*(3*x^2+6*x*y-14*x+3*y^2-14*...
         y+19))*((2*x-3*y)^2*(12*x^2-36*x*y-32*x+27*y^2+48*y+18)+30)-...
         ((x+y+1)^2*(3*x^2+6*x*y-14*x+3*y^2-14*y+19)+1)*((12*x-18*y)*...
         (12*x^2-36*x*y-32*x+27*y^2+48*y+18)-(2*x-3*y)^2*(54*y-36*x+48))];
end

function z=gpf_gfunv(x,y)
    z = [((6.*x+6.*y-14).*(x+y+1).^2+(2.*x+2.*y+2).*(3.*x.^2+6.*x.*y-...
        14.*x+3.*y.^2-14.*y+19)).*((2.*x-3.*y).^2.*(12.*x.^2-36.*x.*y-...
        32.*x+27.*y.^2+48.*y+18)+30)+((x+y+1).^2.*(3.*x.^2+6.*x.*y-14.*...
        x+3.*y.^2-14.*y+19)+1).*((8.*x-12.*y).*(12.*x.^2-36.*x.*y-32.*x+...
        27.*y.^2+48.*y+18)-(2.*x-3.*y).^2.*(36.*y-24.*x+32));
        ((6.*x+6.*y-14).*(x+y+1).^2+(2.*x+2.*y+2).*(3.*x.^2+6.*x.*y-14.*...
        x+3.*y.^2-14.*y+19)).*((2.*x-3.*y).^2.*(12.*x.^2-36.*x.*y-32.*x+...
        27.*y.^2+48.*y+18)+30)-((x+y+1).^2.*(3.*x.^2+6.*x.*y-14.*x+3.*...
        y.^2-14.*y+19)+1).*((12.*x-18.*y).*(12.*x.^2-36.*x.*y-32.*x+27.*...
        y.^2+48.*y+18)-(2.*x-3.*y).^2.*(54.*y-36.*x+48))];
end

function z=gpf_hfun(x,y)
    z=[2*((8*x-12*y)*(12*x^2-36*x*y-32*x+27*y^2+48*y+18)-(2*x-3*y)^2*...
       (36*y-24*x+32))*((6*x+6*y-14)*(x+y+1)^2+(2*x+2*y+2)*(3*x^2+6*x*...
       y-14*x+3*y^2-14*y+19))+((2*x-3*y)^2*(12*x^2-36*x*y-32*x+27*y^2+...
       48*y+18)+30)*(12*x*y-28*y-28*x+2*(2*x+2*y+2)*(6*x+6*y-14)+6*(x+...
       y+1)^2+6*x^2+6*y^2+38)+((x+y+1)^2*(3*x^2+6*x*y-14*x+3*y^2-14*y+...
       19)+1)*(384*y-256*x-288*x*y+24*(2*x-3*y)^2-2*(8*x-12*y)*(36*y-...
       24*x+32)+96*x^2+216*y^2+144),...
       ((8*x-12*y)*(12*x^2-36*x*y-32*x+...
       27*y^2+48*y+18)-(2*x-3*y)^2*(36*y-24*x+32))*((6*x+6*y-14)*(x+y+...
       1)^2+(2*x+2*y+2)*(3*x^2+6*x*y-14*x+3*y^2-14*y+19))-((x+y+1)^2*(3*...
       x^2+6*x*y-14*x+3*y^2-14*y+19)+1)*(576*y-384*x-432*x*y+36*(2*x-3*...
       y)^2-(12*x-18*y)*(36*y-24*x+32)-(8*x-12*y)*(54*y-36*x+48)+144*...
       x^2+324*y^2+216)-((12*x-18*y)*(12*x^2-36*x*y-32*x+27*y^2+48*y+...
       18)-(2*x-3*y)^2*(54*y-36*x+48))*((6*x+6*y-14)*(x+y+1)^2+(2*x+2*y+...
       2)*(3*x^2+6*x*y-14*x+3*y^2-14*y+19))+((2*x-3*y)^2*(12*x^2-36*x*y-...
       32*x+27*y^2+48*y+18)+30)*(12*x*y-28*y-28*x+2*(2*x+2*y+2)*(6*x+6*...
       y-14)+6*(x+y+1)^2+6*x^2+6*y^2+38);
       ((8*x-12*y)*(12*x^2-36*x*y-32*x+27*y^2+48*y+18)-(2*x-3*y)^2*(36*...
       y-24*x+32))*((6*x+6*y-14)*(x+y+1)^2+(2*x+2*y+2)*(3*x^2+6*x*y-14*...
       x+3*y^2-14*y+19))-((x+y+1)^2*(3*x^2+6*x*y-14*x+3*y^2-14*y+19)+...
       1)*(576*y-384*x-432*x*y+36*(2*x-3*y)^2-(12*x-18*y)*(36*y-24*x+...
       32)-(8*x-12*y)*(54*y-36*x+48)+144*x^2+324*y^2+216)-((12*x-18*y)*...
       (12*x^2-36*x*y-32*x+27*y^2+48*y+18)-(2*x-3*y)^2*(54*y-36*x+48))*...
       ((6*x+6*y-14)*(x+y+1)^2+(2*x+2*y+2)*(3*x^2+6*x*y-14*x+3*y^2-14*...
       y+19))+((2*x-3*y)^2*(12*x^2-36*x*y-32*x+27*y^2+48*y+18)+30)*(12*...
       x*y-28*y-28*x+2*(2*x+2*y+2)*(6*x+6*y-14)+6*(x+y+1)^2+6*x^2+6*y^2+...
       38),...
       ((2*x-3*y)^2*(12*x^2-36*x*y-32*x+27*y^2+48*y+18)+30)*(12*x*y-...
       28*y-28*x+2*(2*x+2*y+2)*(6*x+6*y-14)+6*(x+y+1)^2+6*x^2+6*y^2+38)-...
       2*((12*x-18*y)*(12*x^2-36*x*y-32*x+27*y^2+48*y+18)-(2*x-3*y)^2*...
       (54*y-36*x+48))*((6*x+6*y-14)*(x+y+1)^2+(2*x+2*y+2)*(3*x^2+6*x*...
       y-14*x+3*y^2-14*y+19))+((x+y+1)^2*(3*x^2+6*x*y-14*x+3*y^2-14*y+...
       19)+1)*(864*y-576*x-648*x*y+54*(2*x-3*y)^2-2*(12*x-18*y)*(54*y-...
       36*x+48)+216*x^2+486*y^2+324)];
end

function z=gpf_hfunv(x,y)
    z=[2.*((8.*x-12.*y).*(12.*x.^2-36.*x.*y-32.*x+27.*y.^2+48.*y+18)-(2.*x-3.*y).^2.*...
       (36.*y-24.*x+32)).*((6.*x+6.*y-14).*(x+y+1).^2+(2.*x+2.*y+2).*(3.*x.^2+6.*x.*...
       y-14.*x+3.*y.^2-14.*y+19))+((2.*x-3.*y).^2.*(12.*x.^2-36.*x.*y-32.*x+27.*y.^2+...
       48.*y+18)+30).*(12.*x.*y-28.*y-28.*x+2.*(2.*x+2.*y+2).*(6.*x+6.*y-14)+6.*(x+...
       y+1).^2+6.*x.^2+6.*y.^2+38)+((x+y+1).^2.*(3.*x.^2+6.*x.*y-14.*x+3.*y.^2-14.*y+...
       19)+1).*(384.*y-256.*x-288.*x.*y+24.*(2.*x-3.*y).^2-2.*(8.*x-12.*y).*(36.*y-...
       24.*x+32)+96.*x.^2+216.*y.^2+144),...
       ((8.*x-12.*y).*(12.*x.^2-36.*x.*y-32.*x+...
       27.*y.^2+48.*y+18)-(2.*x-3.*y).^2.*(36.*y-24.*x+32)).*((6.*x+6.*y-14).*(x+y+...
       1).^2+(2.*x+2.*y+2).*(3.*x.^2+6.*x.*y-14.*x+3.*y.^2-14.*y+19))-((x+y+1).^2.*(3.*...
       x.^2+6.*x.*y-14.*x+3.*y.^2-14.*y+19)+1).*(576.*y-384.*x-432.*x.*y+36.*(2.*x-3.*...
       y).^2-(12.*x-18.*y).*(36.*y-24.*x+32)-(8.*x-12.*y).*(54.*y-36.*x+48)+144.*...
       x.^2+324.*y.^2+216)-((12.*x-18.*y).*(12.*x.^2-36.*x.*y-32.*x+27.*y.^2+48.*y+...
       18)-(2.*x-3.*y).^2.*(54.*y-36.*x+48)).*((6.*x+6.*y-14).*(x+y+1).^2+(2.*x+2.*y+...
       2).*(3.*x.^2+6.*x.*y-14.*x+3.*y.^2-14.*y+19))+((2.*x-3.*y).^2.*(12.*x.^2-36.*x.*y-...
       32.*x+27.*y.^2+48.*y+18)+30).*(12.*x.*y-28.*y-28.*x+2.*(2.*x+2.*y+2).*(6.*x+6.*...
       y-14)+6.*(x+y+1).^2+6.*x.^2+6.*y.^2+38);
       ((8.*x-12.*y).*(12.*x.^2-36.*x.*y-32.*x+27.*y.^2+48.*y+18)-(2.*x-3.*y).^2.*(36.*...
       y-24.*x+32)).*((6.*x+6.*y-14).*(x+y+1).^2+(2.*x+2.*y+2).*(3.*x.^2+6.*x.*y-14.*...
       x+3.*y.^2-14.*y+19))-((x+y+1).^2.*(3.*x.^2+6.*x.*y-14.*x+3.*y.^2-14.*y+19)+...
       1).*(576.*y-384.*x-432.*x.*y+36.*(2.*x-3.*y).^2-(12.*x-18.*y).*(36.*y-24.*x+...
       32)-(8.*x-12.*y).*(54.*y-36.*x+48)+144.*x.^2+324.*y.^2+216)-((12.*x-18.*y).*...
       (12.*x.^2-36.*x.*y-32.*x+27.*y.^2+48.*y+18)-(2.*x-3.*y).^2.*(54.*y-36.*x+48)).*...
       ((6.*x+6.*y-14).*(x+y+1).^2+(2.*x+2.*y+2).*(3.*x.^2+6.*x.*y-14.*x+3.*y.^2-14.*...
       y+19))+((2.*x-3.*y).^2.*(12.*x.^2-36.*x.*y-32.*x+27.*y.^2+48.*y+18)+30).*(12.*...
       x.*y-28.*y-28.*x+2.*(2.*x+2.*y+2).*(6.*x+6.*y-14)+6.*(x+y+1).^2+6.*x.^2+6.*y.^2+...
       38),...
       ((2.*x-3.*y).^2.*(12.*x.^2-36.*x.*y-32.*x+27.*y.^2+48.*y+18)+30).*(12.*x.*y-...
       28.*y-28.*x+2.*(2.*x+2.*y+2).*(6.*x+6.*y-14)+6.*(x+y+1).^2+6.*x.^2+6.*y.^2+38)-...
       2.*((12.*x-18.*y).*(12.*x.^2-36.*x.*y-32.*x+27.*y.^2+48.*y+18)-(2.*x-3.*y).^2.*...
       (54.*y-36.*x+48)).*((6.*x+6.*y-14).*(x+y+1).^2+(2.*x+2.*y+2).*(3.*x.^2+6.*x.*...
       y-14.*x+3.*y.^2-14.*y+19))+((x+y+1).^2.*(3.*x.^2+6.*x.*y-14.*x+3.*y.^2-14.*y+...
       19)+1).*(864.*y-576.*x-648.*x.*y+54.*(2.*x-3.*y).^2-2.*(12.*x-18.*y).*(54.*y-...
       36.*x+48)+216.*x.^2+486.*y.^2+324)];
    z = reshape_hessian(z);
end


% ACKLEY FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z=ackley_fun(x,y)
    z=-20*exp(-0.2*sqrt(0.5*(x^2+y^2)))-exp(0.5*(cos(2*pi*x)+...
        cos(2*pi*y)))+exp(1)+20;
end

function z=ackley_funv(x,y)
    z=-20.*exp(-0.2.*sqrt(0.5.*(x.^2+y.^2)))-exp(0.5.*(cos(2.*pi.*x)+...
        cos(2.*pi.*y)))+exp(1)+20;
end

function z=ackley_gfun(x,y)
        z=[pi*exp(cos(2*pi*x)/2+cos(2*pi*y)/2)*sin(2*pi*x)+...
            (2*x*exp(-(x^2/2+y^2/2)^(1/2)/5))/(x^2/2+y^2/2)^(1/2);
           pi*exp(cos(2*pi*x)/2+cos(2*pi*y)/2)*sin(2*pi*y)+...
           (2*y*exp(-(x^2/2+y^2/2)^(1/2)/5))/(x^2/2+y^2/2)^(1/2)];
end

function z=ackley_gfunv(x,y)
        z=[pi.*exp(cos(2.*pi.*x)./2+cos(2.*pi.*y)./2).*sin(2.*pi.*x)+...
           (2.*x.*exp(-(x.^2./2+y.^2./2).^(1./2)./5))./(x.^2./2+y.^2./2).^(1./2);
           pi.*exp(cos(2.*pi.*x)./2+cos(2.*pi.*y)./2).*sin(2.*pi.*y)+...
           (2.*y.*exp(-(x.^2./2+y.^2./2).^(1./2)./5))./(x.^2./2+y.^2./2).^(1./2)]; 
end

function z=ackley_hfun(x,y)
    if x==0.0&&y==0.0
        %this isnt the correct answer, but it will give indicate that the
        %solution is positive definite
        z=[1,0;0,1]; 
    else
    z=[(2*exp(-(x^2/2+y^2/2)^(1/2)/5))/(x^2/2+y^2/2)^(1/2)+2*pi^2*...
        exp(cos(2*pi*x)/2+cos(2*pi*y)/2)*cos(2*pi*x)-(x^2*exp(-(x^2/2+...
        y^2/2)^(1/2)/5))/(5*(x^2/2+y^2/2))-(x^2*exp(-(x^2/2+...
        y^2/2)^(1/2)/5))/(x^2/2+y^2/2)^(3/2)-pi^2*exp(cos(2*pi*x)/2+...
        cos(2*pi*y)/2)*sin(2*pi*x)^2,-(x*y*exp(-(x^2/2+y^2/2)^(1/2)/...
        5))/(5*(x^2/2+y^2/2))-(x*y*exp(-(x^2/2+y^2/2)^(1/2)/5))/(x^2/2+...
        y^2/2)^(3/2)-pi^2*exp(cos(2*pi*x)/2+cos(2*pi*y)/2)*...
        sin(2*pi*x)*sin(2*pi*y);
       -(x*y*exp(-(x^2/2+y^2/2)^(1/2)/5))/(5*(x^2/2+y^2/2))-(x*y*...
       exp(-(x^2/2+y^2/2)^(1/2)/5))/(x^2/2+y^2/2)^(3/2)-pi^2*...
       exp(cos(2*pi*x)/2+cos(2*pi*y)/2)*sin(2*pi*x)*sin(2*pi*y),(2*...
       exp(-(x^2/2+y^2/2)^(1/2)/5))/(x^2/2+y^2/2)^(1/2)+2*pi^2*...
       exp(cos(2*pi*x)/2+cos(2*pi*y)/2)*cos(2*pi*y)-(y^2*exp(-(x^2/2+...
       y^2/2)^(1/2)/5))/(5*(x^2/2+y^2/2))-(y^2*exp(-(x^2/2+y^2/2)^(1/2)/...
       5))/(x^2/2+y^2/2)^(3/2)-pi^2*exp(cos(2*pi*x)/2+...
       cos(2*pi*y)/2)*sin(2*pi*y)^2];
    end
end

function z=ackley_hfunv(x,y)
    z=[(2.*exp(-(x.^2./2+y.^2./2).^(1./2)./5))./(x.^2./2+y.^2./2).^(1./2)+2.*pi.^2.*...
       exp(cos(2.*pi.*x)./2+cos(2.*pi.*y)./2).*cos(2.*pi.*x)-(x.^2.*exp(-(x.^2./2+...
       y.^2./2).^(1./2)./5))./(5.*(x.^2./2+y.^2./2))-(x.^2.*exp(-(x.^2./2+...
       y.^2./2).^(1./2)./5))./(x.^2./2+y.^2./2).^(3./2)-pi.^2.*exp(cos(2.*pi.*x)./2+...
       cos(2.*pi.*y)./2).*sin(2.*pi.*x).^2,-(x.*y.*exp(-(x.^2./2+y.^2./2).^(1./2)./...
       5))./(5.*(x.^2./2+y.^2./2))-(x.*y.*exp(-(x.^2./2+y.^2./2).^(1./2)./5))./(x.^2./2+...
       y.^2./2).^(3./2)-pi.^2.*exp(cos(2.*pi.*x)./2+cos(2.*pi.*y)./2).*...
       sin(2.*pi.*x).*sin(2.*pi.*y);
       -(x.*y.*exp(-(x.^2./2+y.^2./2).^(1./2)./5))./(5.*(x.^2./2+y.^2./2))-(x.*y.*...
       exp(-(x.^2./2+y.^2./2).^(1./2)./5))./(x.^2./2+y.^2./2).^(3./2)-pi.^2.*...
       exp(cos(2.*pi.*x)./2+cos(2.*pi.*y)./2).*sin(2.*pi.*x).*sin(2.*pi.*y),(2.*...
       exp(-(x.^2./2+y.^2./2).^(1./2)./5))./(x.^2./2+y.^2./2).^(1./2)+2.*pi.^2.*...
       exp(cos(2.*pi.*x)./2+cos(2.*pi.*y)./2).*cos(2.*pi.*y)-(y.^2.*exp(-(x.^2./2+...
       y.^2./2).^(1./2)./5))./(5.*(x.^2./2+y.^2./2))-(y.^2.*exp(-(x.^2./2+y.^2./2).^(1./2)./...
       5))./(x.^2./2+y.^2./2).^(3./2)-pi.^2.*exp(cos(2.*pi.*x)./2+...
       cos(2.*pi.*y)./2).*sin(2.*pi.*y).^2];
    z = reshape_hessian(z);

end


% RASTRIGIN FUNCTION 4 DIMENSIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z = rastrigin4_fun(w1,w2,w3,w4)
    z = w1^2-10*cos(2*pi*w2)-10*cos(2*pi*w3)-10*cos(2*pi*w4)-...
        10*cos(2*pi*w1)+w2^2+w3^2+w4^2+40;
end

function z = rastrigin4_funv(w1,w2,w3,w4)
    z = w1.^2-10.*cos(2.*pi.*w2)-10.*cos(2.*pi.*w3)-10.*cos(2.*pi.*w4)-...
        10.*cos(2.*pi.*w1)+w2.^2+w3.^2+w4.^2+40;
end

function z = rastrigin4_gfun(w1,w2,w3,w4)
    z = [2*w1 + 20*pi*sin(2*pi*w1);
         2*w2 + 20*pi*sin(2*pi*w2);
         2*w3 + 20*pi*sin(2*pi*w3);
         2*w4 + 20*pi*sin(2*pi*w4)];
end

function z = rastrigin4_gfunv(w1,w2,w3,w4)
    z = [2.*w1 + 20.*pi.*sin(2.*pi.*w1);
         2.*w2 + 20.*pi.*sin(2.*pi.*w2);
         2.*w3 + 20.*pi.*sin(2.*pi.*w3);
         2.*w4 + 20.*pi.*sin(2.*pi.*w4)];
end

function z = rastrigin4_hfun(w1,w2,w3,w4)
    z = [40*pi^2*cos(2*pi*w1)+2,0,0,0;
         0,40*pi^2*cos(2*pi*w2)+2,0,0;
         0,0,40*pi^2*cos(2*pi*w3)+2,0;
         0,0,0,40*pi^2*cos(2*pi*w4)+2];
end

function z = rastrigin4_hfunv(w1,w2,w3,w4)
    z = [40.*pi.^2.*cos(2.*pi.*w1)+2,0,0,0;
         0,40.*pi.^2.*cos(2.*pi.*w2)+2,0,0;
         0,0,40.*pi.^2.*cos(2.*pi.*w3)+2,0;
         0,0,0,40.*pi.^2.*cos(2.*pi.*w4)+2];
    z = reshape_hessian(z);
end

% RASTRIGIN FUNCTION 2 DIMENSIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z = rastrigin2_fun(w1,w2)
    z = w1^2-10*cos(2*pi*w2)-10*cos(2*pi*w1)+w2^2+20;
end

function z = rastrigin2_funv(w1,w2)
    z = w1.^2-10.*cos(2.*pi.*w2)-10.*cos(2.*pi.*w1)+w2.^2+20;
end

function z = rastrigin2_gfun(w1,w2)
    z = [2*w1 + 20*pi*sin(2*pi*w1);
         2*w2 + 20*pi*sin(2*pi*w2)];
end

function z = rastrigin2_gfunv(w1,w2)
    z = [2.*w1 + 20.*pi.*sin(2.*pi.*w1);
         2.*w2 + 20.*pi.*sin(2.*pi.*w2)];
end

function z = rastrigin2_hfun(w1,w2)
    z = [40*pi^2*cos(2*pi*w1)+2,0;
         0,40*pi^2*cos(2*pi*w2)+2];
end

function z = rastrigin2_hfunv(w1,w2)
    z = [40.*pi.^2.*cos(2.*pi.*w1)+2,0;
         0,40.*pi.^2.*cos(2.*pi.*w2)+2];
    z = reshape_hessian(z);
end

% SPHERE FUNCTION IN 2 DIMENSIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z = sphere2_fun(w1,w2)
    z = w1^2+w2^2;
end

function z = sphere2_funv(w1,w2)
    z = w1.^2+w2.^2;
end

function z = sphere2_gfun(w1,w2)
    z = [2*w1;
         2*w2];
end

function z = sphere2_gfunv(w1,w2)
    z = [2.*w1;
         2.*w2];
end

function z = sphere2_hfun(~,~)
    z = [2,0;
         0,2];
end

function z = sphere2_hfunv(~,~)
    z = [2,0;
         0,2];
end

% SPHERE FUNCTION IN 4 DIMENSIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z = sphere4_fun(w1,w2,w3,w4)
    z = w1^2+w2^2+w3^2+w4^2;
end

function z = sphere4_funv(w1,w2,w3,~)
    z = w1.^2+w2.^2+w3.^2+w3.^2;
end

function z = sphere4_gfun(w1,w2,w3,w4)
    z = [2*w1;
         2*w2;
         2*w3;
         2*w4];
end

function z = sphere4_gfunv(w1,w2,w3,w4)
    z = [2.*w1;
         2.*w2;
         2.*w3;
         2.*w4];
end

function z = sphere4_hfun(~,~,~,~)
    z = [2,0,0,0;
         0,2,0,0;
         0,0,2,0;
         0,0,0,2];
end

function z = sphere4_hfunv(~,~,~,~)
    z = [2,0,0,0;
         0,2,0,0;
         0,0,2,0;
         0,0,0,2];
end

% ROSENBROCK FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z = rosenbrock_fun(x,y)
    z = 100*y+(x-1)^2-100*x^2;
end

function z = rosenbrock_funv(x,y)
    z = 100.*y+(x-1).^2-100.*x.^2;
end

function z = rosenbrock_gfun(x,~)
    z = [-198*x-2;
        100];
end

function z = rosenbrock_gfunv(x,~)
    z = [-198.*x-2;
        100];
end

function z = rosenbrock_hfun(~,~)
    z = [-198,0;
        0,0];
end

function z = rosenbrock_hfunv(~,~)
    z = [-198,0;
        0,0];
end

% BOOTH FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z = booth_fun(x,y)
    z = (2*x+y-5)^2+(x+2*y-7)^2;
end

function z = booth_funv(x,y)
    z = (2.*x+y-5).^2+(x+2.*y-7).^2;
end

function z = booth_gfun(x,y)
    z = [10*x+8*y-34;
          8*x+10*y-38];
end

function z = booth_gfunv(x,y)
    z = [10.*x+8.*y-34;
          8.*x+10.*y-38];
end

function z = booth_hfun(~,~)
    z = [10,8;
        8,10];
end

function z = booth_hfunv(~,~)
    z = [10,8;
        8,10];
end

% BULKIN No. 6 FUCNTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z = bulkinN6_fun(x,y)
    z = abs(x+10)/100+100*abs(x^2/100-y)^(1/2);
end

function z = bulkinN6_funv(x,y)
    z = abs(x+10)./100+100.*abs(x.^2./100-y).^(1./2);
end

function z = bulkinN6_gfun(x,y)
    z = [sign(x+10)/100+(x*sign(x^2/100-y))/abs(x^2/100-y)^(1/2);
                 -(50*sign(x^2/100-y))/abs(x^2/100-y)^(1/2)];
end

function z = bulkinN6_gfunv(x,y)
    z = [sign(x+10)./100+(x.*sign(x.^2./100-y))./abs(x.^2./100-y).^(1./2);
                 -(50.*sign(x.^2./100-y))./abs(x.^2./100-y).^(1./2)];
end

function z = bulkinN6_hfun(x,y)
    z = [dirac(x+10)/50+sign(x^2/100-y)/abs(x^2/100-y)^(1/2)+(x^2*...
        dirac(x^2/100-y))/(25*abs(x^2/100-y)^(1/2))-(x^2*sign(x^2/100-...
        y)^2)/(100*abs(x^2/100-y)^(3/2)),...
        (x*sign(x^2/100-y)^2)/(2*abs(x^2/100-y)^(3/2))-(2*x*...
        dirac(x^2/100-y))/abs(x^2/100-y)^(1/2);
         (x*sign(x^2/100-y)^2)/(2*abs(x^2/100-y)^(3/2))-(2*x*...
         dirac(x^2/100-y))/abs(x^2/100-y)^(1/2),...
         (100*dirac(x^2/100-y))/abs(x^2/100-y)^(1/2)-(25*sign(x^2/100-...
         y)^2)/abs(x^2/100-y)^(3/2)];
end

function z = bulkinN6_hfunv(x,y)
    z = [dirac(x+10)./50+sign(x.^2./100-y)./abs(x.^2./100-y).^(1./2)+(x.^2.*...
        dirac(x.^2./100-y))./(25.*abs(x.^2./100-y).^(1./2))-(x.^2.*sign(x.^2./100-...
        y).^2)./(100.*abs(x.^2./100-y).^(3./2)),...
        (x.*sign(x.^2./100-y).^2)./(2.*abs(x.^2./100-y).^(3./2))-(2.*x.*...
        dirac(x.^2./100-y))./abs(x.^2./100-y).^(1./2);
         (x.*sign(x.^2./100-y).^2)./(2.*abs(x.^2./100-y).^(3./2))-(2.*x.*...
         dirac(x.^2./100-y))./abs(x.^2./100-y).^(1./2),...
         (100.*dirac(x.^2./100-y))./abs(x.^2./100-y).^(1./2)-(25.*sign(x.^2./100-...
         y).^2)./abs(x.^2./100-y).^(3./2)];
    z = reshape_hessian(z);
end

% HIMMELBLAU FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z = himmelblau_fun(x,y)
    z = (y^2+x-7)^2+(x^2+y-11)^2;
end

function z = himmelblau_funv(x,y)
    z = (y.^2+x-7).^2+(x.^2+y-11).^2;
end

function z = himmelblau_gfun(x,y)
    z = [2*x+4*x*(x^2+y-11)+2*y^2-14;
         2*y+4*y*(y^2+x-7)+2*x^2-22];
end

function z = himmelblau_gfunv(x,y)
    z = [2.*x+4.*x.*(x.^2+y-11)+2.*y.^2-14;
         2.*y+4.*y.*(y.^2+x-7)+2.*x.^2-22];
end

function z = himmelblau_hfun(x,y)
    z = [12*x^2+4*y-42,4*x+4*y;
        4*x+4*y,12*y^2+4*x-26];
end

function z = himmelblau_hfunv(x,y)
    z = [12.*x.^2+4.*y-42,4.*x+4.*y;
        4.*x+4.*y,12.*y.^2+4.*x-26];
    z = reshape_hessian(z);
end

% HOLDER TABLE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function z = holdertable_fun(x,y)
    z = -exp(abs((x^2+y^2)^(1/2)/pi-1))*abs(cos(y)*sin(x));
end

function z = holdertable_funv(x,y)
    z = -exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*abs(cos(y).*sin(x));
end

function z = holdertable_gfun(x,y)
    z = [-exp(abs((x^2+y^2)^(1/2)/pi-1))*sign(cos(y)*sin(x))*cos(x)*...
         cos(y)-(x*sign((x^2+y^2)^(1/2)/pi-1)*exp(abs((x^2+y^2)^(1/2)/...
         pi-1))*abs(cos(y)*sin(x)))/(pi*(x^2+y^2)^(1/2));
         exp(abs((x^2+y^2)^(1/2)/pi-1))*sign(cos(y)*sin(x))*sin(x)*...
         sin(y)-(y*sign((x^2+y^2)^(1/2)/pi-1)*exp(abs((x^2+y^2)^(1/2)/...
         pi-1))*abs(cos(y)*sin(x)))/(pi*(x^2+y^2)^(1/2))];
end

function z = holdertable_gfunv(x,y)
    z = [-exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*sign(cos(y).*sin(x)).*cos(x).*...
         cos(y)-(x.*sign((x.^2+y.^2).^(1./2)./pi-1).*exp(abs((x.^2+y.^2).^(1./2)./...
         pi-1)).*abs(cos(y).*sin(x)))./(pi.*(x.^2+y.^2).^(1./2));
         exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*sign(cos(y).*sin(x)).*sin(x).*...
         sin(y)-(y.*sign((x.^2+y.^2).^(1./2)./pi-1).*exp(abs((x.^2+y.^2).^(1./2)./...
         pi-1)).*abs(cos(y).*sin(x)))./(pi.*(x.^2+y.^2).^(1./2))];
end

function z = holdertable_hfun(x,y)
    z = [exp(abs((x^2+y^2)^(1/2)/pi-1))*sign(cos(y)*sin(x))*cos(y)*...
        sin(x)-2*exp(abs((x^2+y^2)^(1/2)/pi-1))*dirac(cos(y)*sin(x))*...
        cos(x)^2*cos(y)^2-(sign((x^2+y^2)^(1/2)/pi-1)*exp(abs((x^2+y^2)^...
        (1/2)/pi-1))*abs(cos(y)*sin(x)))/(pi*(x^2+y^2)^(1/2))-(x^2*...
        sign((x^2+y^2)^(1/2)/pi-1)^2*exp(abs((x^2+y^2)^(1/2)/pi-1))*...
        abs(cos(y)*sin(x)))/(pi^2*(x^2+y^2))-(2*x^2*dirac((x^2+y^2)^...
        (1/2)/pi-1)*exp(abs((x^2+y^2)^(1/2)/pi-1))*abs(cos(y)*sin(x)))/...
        (pi^2*(x^2+y^2))+(x^2*sign((x^2+y^2)^(1/2)/pi-1)*exp(abs((x^2+...
        y^2)^(1/2)/pi-1))*abs(cos(y)*sin(x)))/(pi*(x^2+y^2)^(3/2))-(2*x*...
        sign((x^2+y^2)^(1/2)/pi-1)*exp(abs((x^2+y^2)^(1/2)/pi-1))*...
        sign(cos(y)*sin(x))*cos(x)*cos(y))/(pi*(x^2+y^2)^(1/2)),...
        exp(abs((x^2+y^2)^(1/2)/pi-1))*sign(cos(y)*sin(x))*cos(x)*...
        sin(y)+2*exp(abs((x^2+y^2)^(1/2)/pi-1))*dirac(cos(y)*sin(x))*...
        cos(x)*cos(y)*sin(x)*sin(y)-(x*y*sign((x^2+y^2)^(1/2)/pi-1)^2*...
        exp(abs((x^2+y^2)^(1/2)/pi-1))*abs(cos(y)*sin(x)))/(pi^2*(x^2+...
        y^2))-(2*x*y*dirac((x^2+y^2)^(1/2)/pi-1)*exp(abs((x^2+y^2)^...
        (1/2)/pi-1))*abs(cos(y)*sin(x)))/(pi^2*(x^2+y^2))+(x*y*sign((x^...
        2+y^2)^(1/2)/pi-1)*exp(abs((x^2+y^2)^(1/2)/pi-1))*abs(cos(y)*...
        sin(x)))/(pi*(x^2+y^2)^(3/2))-(y*sign((x^2+y^2)^(1/2)/pi-1)*...
        exp(abs((x^2+y^2)^(1/2)/pi-1))*sign(cos(y)*sin(x))*cos(x)*...
        cos(y))/(pi*(x^2+y^2)^(1/2))+(x*sign((x^2+y^2)^(1/2)/pi-1)*...
        exp(abs((x^2+y^2)^(1/2)/pi-1))*sign(cos(y)*sin(x))*sin(x)*...
        sin(y))/(pi*(x^2+y^2)^(1/2));
        exp(abs((x^2+y^2)^(1/2)/pi-1))*sign(cos(y)*sin(x))*cos(x)*...
        sin(y)+2*exp(abs((x^2+y^2)^(1/2)/pi-1))*dirac(cos(y)*sin(x))*...
        cos(x)*cos(y)*sin(x)*sin(y)-(x*y*sign((x^2+y^2)^(1/2)/pi-1)^2*...
        exp(abs((x^2+y^2)^(1/2)/pi-1))*abs(cos(y)*sin(x)))/(pi^2*(x^2+...
        y^2))-(2*x*y*dirac((x^2+y^2)^(1/2)/pi-1)*exp(abs((x^2+y^2)^...
        (1/2)/pi-1))*abs(cos(y)*sin(x)))/(pi^2*(x^2+y^2))+(x*y*sign((x^...
        2+y^2)^(1/2)/pi-1)*exp(abs((x^2+y^2)^(1/2)/pi-1))*abs(cos(y)*...
        sin(x)))/(pi*(x^2+y^2)^(3/2))-(y*sign((x^2+y^2)^(1/2)/pi-1)*...
        exp(abs((x^2+y^2)^(1/2)/pi-1))*sign(cos(y)*sin(x))*cos(x)*...
        cos(y))/(pi*(x^2+y^2)^(1/2))+(x*sign((x^2+y^2)^(1/2)/pi-1)*...
        exp(abs((x^2+y^2)^(1/2)/pi-1))*sign(cos(y)*sin(x))*sin(x)*...
        sin(y))/(pi*(x^2+y^2)^(1/2)),...
        exp(abs((x^2+y^2)^(1/2)/pi-1))*sign(cos(y)*sin(x))*cos(y)*...
        sin(x)-2*exp(abs((x^2+y^2)^(1/2)/pi-1))*dirac(cos(y)*sin(x))*...
        sin(x)^2*sin(y)^2-(sign((x^2+y^2)^(1/2)/pi-1)*exp(abs((x^2+...
        y^2)^(1/2)/pi-1))*abs(cos(y)*sin(x)))/(pi*(x^2+y^2)^(1/2))-(y^2*...
        sign((x^2+y^2)^(1/2)/pi-1)^2*exp(abs((x^2+y^2)^(1/2)/pi-1))*...
        abs(cos(y)*sin(x)))/(pi^2*(x^2+y^2))-(2*y^2*dirac((x^2+y^2)^...
        (1/2)/pi-1)*exp(abs((x^2+y^2)^(1/2)/pi-1))*abs(cos(y)*sin(x)))/...
        (pi^2*(x^2+y^2))+(y^2*sign((x^2+y^2)^(1/2)/pi-1)*exp(abs((x^2+...
        y^2)^(1/2)/pi-1))*abs(cos(y)*sin(x)))/(pi*(x^2+y^2)^(3/2))+(2*y*...
        sign((x^2+y^2)^(1/2)/pi-1)*exp(abs((x^2+y^2)^(1/2)/pi-1))*...
        sign(cos(y)*sin(x))*sin(x)*sin(y))/(pi*(x^2+y^2)^(1/2))];
end

function z = holdertable_hfunv(x,y)
    z = [exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*sign(cos(y).*sin(x)).*...
        cos(y).*sin(x)-2.*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*...
        dirac(cos(y).*sin(x)).*cos(x).^2.*cos(y).^2-(sign((x.^2+y.^2).^...
        (1./2)./pi-1).*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*abs(cos(y).*...
        sin(x)))./(pi.*(x.^2+y.^2).^(1./2))-(x.^2.*sign((x.^2+y.^2).^...
        (1./2)./pi-1).^2.*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*...
        abs(cos(y).*sin(x)))./(pi.^2.*(x.^2+y.^2))-(2.*x.^2.*dirac((x.^...
        2+y.^2).^(1./2)./pi-1).*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*...
        abs(cos(y).*sin(x)))./(pi.^2.*(x.^2+y.^2))+(x.^2.*sign((x.^2+y.^...
        2).^(1./2)./pi-1).*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*...
        abs(cos(y).*sin(x)))./(pi.*(x.^2+y.^2).^(3./2))-(2.*x.*sign((x.^...
        2+y.^2).^(1./2)./pi-1).*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*...
sign(cos(y).*sin(x)).*cos(x).*cos(y))./(pi.*(x.^2+y.^2).^(1./2)),...
exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*sign(cos(y).*sin(x)).*cos(x).*...
sin(y)+2.*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*dirac(cos(y).*sin(x)).*...
cos(x).*cos(y).*sin(x).*sin(y)-(x.*y.*sign((x.^2+y.^2).^(1./2)./pi-1).^2.*...
exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*abs(cos(y).*sin(x)))./(pi.^2.*(x.^2+...
y.^2))-(2.*x.*y.*dirac((x.^2+y.^2).^(1./2)./pi-1).*exp(abs((x.^2+y.^2).^...
(1./2)./pi-1)).*abs(cos(y).*sin(x)))./(pi.^2.*(x.^2+y.^2))+(x.*y.*sign((x.^...
2+y.^2).^(1./2)./pi-1).*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*abs(cos(y).*...
sin(x)))./(pi.*(x.^2+y.^2).^(3./2))-(y.*sign((x.^2+y.^2).^(1./2)./pi-1).*...
exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*sign(cos(y).*sin(x)).*cos(x).*...
cos(y))./(pi.*(x.^2+y.^2).^(1./2))+(x.*sign((x.^2+y.^2).^(1./2)./pi-1).*...
exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*sign(cos(y).*sin(x)).*sin(x).*...
sin(y))./(pi.*(x.^2+y.^2).^(1./2));
exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*sign(cos(y).*sin(x)).*cos(x).*...
sin(y)+2.*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*dirac(cos(y).*sin(x)).*...
cos(x).*cos(y).*sin(x).*sin(y)-(x.*y.*sign((x.^2+y.^2).^(1./2)./pi-1).^2.*...
exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*abs(cos(y).*sin(x)))./(pi.^2.*(x.^2+...
y.^2))-(2.*x.*y.*dirac((x.^2+y.^2).^(1./2)./pi-1).*exp(abs((x.^2+y.^2).^...
(1./2)./pi-1)).*abs(cos(y).*sin(x)))./(pi.^2.*(x.^2+y.^2))+(x.*y.*sign((x.^...
2+y.^2).^(1./2)./pi-1).*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*abs(cos(y).*...
sin(x)))./(pi.*(x.^2+y.^2).^(3./2))-(y.*sign((x.^2+y.^2).^(1./2)./pi-1).*...
exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*sign(cos(y).*sin(x)).*cos(x).*...
cos(y))./(pi.*(x.^2+y.^2).^(1./2))+(x.*sign((x.^2+y.^2).^(1./2)./pi-1).*...
exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*sign(cos(y).*sin(x)).*sin(x).*...
sin(y))./(pi.*(x.^2+y.^2).^(1./2)),...
exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*sign(cos(y).*sin(x)).*cos(y).*...
sin(x)-2.*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*dirac(cos(y).*sin(x)).*...
sin(x).^2.*sin(y).^2-(sign((x.^2+y.^2).^(1./2)./pi-1).*exp(abs((x.^2+...
y.^2).^(1./2)./pi-1)).*abs(cos(y).*sin(x)))./(pi.*(x.^2+y.^2).^(1./2))-(y.^2.*...
sign((x.^2+y.^2).^(1./2)./pi-1).^2.*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*...
abs(cos(y).*sin(x)))./(pi.^2.*(x.^2+y.^2))-(2.*y.^2.*dirac((x.^2+y.^2).^...
(1./2)./pi-1).*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*abs(cos(y).*sin(x)))./...
(pi.^2.*(x.^2+y.^2))+(y.^2.*sign((x.^2+y.^2).^(1./2)./pi-1).*exp(abs((x.^2+...
y.^2).^(1./2)./pi-1)).*abs(cos(y).*sin(x)))./(pi.*(x.^2+y.^2).^(3./2))+(2.*y.*...
sign((x.^2+y.^2).^(1./2)./pi-1).*exp(abs((x.^2+y.^2).^(1./2)./pi-1)).*...
sign(cos(y).*sin(x)).*sin(x).*sin(y))./(pi.*(x.^2+y.^2).^(1./2))];
    z = reshape_hessian(z);
end

function mout = reshape_hessian(m)
    l    = length(m)/2;
    mout = zeros(2,2,l);

    for i = 1:l
        mout(:,:,i) = m(1:2,(2*i-1):2*i);
    end

end