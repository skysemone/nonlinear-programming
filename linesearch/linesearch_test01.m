% Trust region minimization test script
format long

%function, gradient, Hessian definitions
f  =@(x) x(1).^4+3*x(1).^3-2*x(1).^2+x(2).^4+9*x(2).^2+x(1).*x(2)+x(3).^2;
df =@(x)[4*x(1).^3+9*x(1).^2-4*x(1)+x(2);4*x(2).^3+18*x(2)+x(1);2*x(3)];
hf =@(x)[12*x(1).^2+18*x(1)-4,1,0;1,12*x(2).^2+18,0;0,0,2];

%initial x point
x0 =[-2.9;-0.80;4];

%run minimization with print-outs enabled
[x_min,x_list]=linesearch_armijo_min(f,df,hf,x0);

%or alternatively for the wolfe criteria...
%[x_min,x_list]=linesearch_wolfe_min(f,df,hf,x0);

%plot the function and the approximations from the trust-region min
%   set x1, x2 values and make a grid
x1 = -3:0.01:-1.0;
x2 = -1.0:0.01:2.0;
[X1,X2] = ndgrid(x1,x2);

%set x3 to be zeros to make it easy to view in a simple contour plot
X3 = zeros(size(X1));

%flip the matricies to be used by function f
v = num2cell([X1(:), X2(:),X3(:)], 2);
Z = reshape(cellfun(f, v), size(X1));

%plot the data in a contour plot
hold on;
contourf(X1,X2,Z,20)
colormap('summer')
plot(x_list(:,1),x_list(:,2),'-k.')
hold off;

