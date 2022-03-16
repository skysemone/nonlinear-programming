%
%
%
%
%
% Symbolic column vector of length 2
y=sym('y',[2,1]);

% Symbolic column vector of length s, 4 is used for many of the benchmark
% functions
s = 4;
w=sym('w',[s,1]);

% Symbolic formulation of the Beale Function returning a double for the
% function, gradient, and hessian.
ff1 = (1.5-y(1)+y(1)*y(2))^2+(2.25-y(1)+y(1)*y(2)^2)^2 +...
    (2.625-y(1)+y(1)*y(2)^3)^2;
f1  = @(y0) double(subs(ff1,y,y0));
g1  = @(y0) double(subs(gradient(ff1,y),y,y0));
h1  = @(y0) double(subs(hessian(ff1,y),y,y0));

% Symbolic formulation of the Goldstein-Price Function returning a double
% for the function, gradient, and hessian.
ff2 = (1+(y(1)+y(2)+1)^2*(19-14*y(1)+3*y(1)^2-14*y(2)+6*y(1)*y(2)+...
    3*y(2)^2))*(30+(2*y(1)-3*y(2))^2*(18-32*y(1)+12*y(1)^2+48*y(2)-...
    36*y(1)*y(2)+27*y(2)^2));

f2  = @(y0) double(subs(ff2,y,y0));
g2  = @(y0) double(subs(gradient(ff2,y),y,y0));
h2  = @(y0) double(subs(hessian(ff2,y),y,y0));

% Symbolic formulation of the Ackley Function returning a double for the
% function, gradient, and hessian.
ff3 = -20*exp(-0.2*sqrt(0.5*(y(1)^2+y(2)^2)))-exp(0.5*(cos(2*pi*y(1))+...
    cos(2*pi*y(2))))+exp(1)+20;

f3  = @(y0) (subs(ff3,y,y0));
g3  = @(y0) (subs(gradient(ff3,y),y,y0));
h3  = @(y0) (subs(hessian(ff3,y),y,y0));

% Symbolic formulation of the Rastrigin Function returning a double for the
% function, gradient, and hessian in s dimensions.  Default dimensions are
% set to 4.
ff4 = 10*s+sum(w.^2-10*cos(2*pi*w));

f4 =@(y0) double(subs(ff4,w,y0));
g4 =@(y0) double(subs(gradient(ff4,w),w,y0));
h4 =@(y0) double(subs(hessian(ff4,w),w,y0));

% Symbolic formulation of the Rosenbrock Function returning a double
% for the function, gradient, and hessian.
ff5 = 100*(y(2)-y(1)^2)+(1-y(1))^2;

f5 =@(y0) (subs(ff5,y,y0));
g5 =@(y0) (subs(gradient(ff5,y),y,y0));
h5 =@(y0) (subs(hessian(ff5,y),y,y0));

% Symbolic formulation of the Booth Function returning a double
% for the function, gradient, and hessian.
ff6 = (y(1)+2*y(2)-7)^2+(2*y(1)+y(2)-5)^2;

f6 =@(y0) (subs(ff6,y,y0));
g6 =@(y0) (subs(gradient(ff6,y),y,y0));
h6 =@(y0) (subs(hessian(ff6,y),y,y0));

% Symbolic formulation of the Himmelblau Function returning a double
% for the function, gradient, and hessian.
ff7 = (y(1)^2+y(2)-11)^2+(y(1)+y(2)^2-7)^2;

f7 =@(y0) (subs(ff7,y,y0));
g7 =@(y0) (subs(gradient(ff7,y),y,y0));
h7 =@(y0) (subs(hessian(ff7,y),y,y0));

% Symbolic formulation of the Bulkin No. 6 Function returning a double
% for the function, gradient, and hessian.
ff8 = 100*sqrt(abs(y(2)-0.01*y(1)^2))+0.01*abs(y(1)+10);

f8 =@(y0) (subs(ff8,y,y0));
g8 =@(y0) (subs(gradient(ff8,y),y,y0));
h8 =@(y0) (subs(hessian(ff8,y),y,y0));

% Symbolic formulation of the Bulkin No. 6 Function returning a double
% for the function, gradient, and hessian.
ff9 = -abs(sin(y(1))*cos(y(2))*exp(abs(1-sqrt(y(1)^2+y(2)^2)/pi)));

f9 =@(y0) (subs(ff9,y,y0));
g9 =@(y0) (subs(gradient(ff9,y),y,y0));
h9 =@(y0) (subs(hessian(ff9,y),y,y0));


% Reset the x and y to non-vector symbols.
syms x y;





